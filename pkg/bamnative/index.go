package bamnative

import (
	"bufio"
	"encoding/binary"
	"fmt"
	"io"
	"math"
	"os"
	"path/filepath"
	"sort"
)

// BAI magic number
var baiMagic = [4]byte{'B', 'A', 'I', 1}

// metaBin is the special BAI pseudo-bin that stores per-reference statistics.
const (
	metaBin              = uint32(37450)
	maximumBAICoordinate = int64(1 << 29)
)

// chunk represents a BAI index chunk (a contiguous range of virtual offsets).
type chunk struct {
	start int64
	end   int64
}

// perRefIndex holds the bin index, linear index, and statistics for one reference.
type perRefIndex struct {
	bins      map[uint32][]chunk
	linear    []int64 // minimum voff for each 16 kb window; 0 means unset
	nMapped   uint64
	nUnmapped uint64
	minVoff   int64 // voffBeg of first record seen for this reference
	maxVoff   int64 // voffEnd of last record seen for this reference
	hasAny    bool  // true if at least one record belongs to this reference
}

// alignedLen returns the number of reference bases consumed by a record's CIGAR
// (operations M, D, N, =, X).
func alignedLen(record *Record) int64 {
	var referenceLength int64
	for _, operation := range record.Cigar {
		switch operation.Op {
		case CigarMatch, CigarDeletion, CigarSkip, CigarEqual, CigarMismatch:
			referenceLength += int64(operation.Len)
		}
	}
	return referenceLength
}

func baiReferenceInterval(record *Record) (int, int, error) {
	begin := int64(record.Pos)
	if begin < 0 || begin >= maximumBAICoordinate {
		return 0, 0, fmt.Errorf("position %d is outside the BAI coordinate range", record.Pos)
	}

	referenceLength := alignedLen(record)
	if referenceLength == 0 {
		referenceLength = 1
	}
	end := begin + referenceLength
	if end <= begin || end > maximumBAICoordinate {
		return 0, 0, fmt.Errorf("alignment end %d is outside the BAI coordinate range", end)
	}
	return int(begin), int(end), nil
}

// reg2bin returns the smallest BAI bin that fully contains [beg, end)
// following the formula in the SAM specification §5.3.
func reg2bin(beg, end int) uint32 {
	end--
	if beg>>14 == end>>14 {
		return uint32(((1<<15)-1)/7 + (beg >> 14))
	}
	if beg>>17 == end>>17 {
		return uint32(((1<<12)-1)/7 + (beg >> 17))
	}
	if beg>>20 == end>>20 {
		return uint32(((1<<9)-1)/7 + (beg >> 20))
	}
	if beg>>23 == end>>23 {
		return uint32(((1<<6)-1)/7 + (beg >> 23))
	}
	if beg>>26 == end>>26 {
		return uint32(((1<<3)-1)/7 + (beg >> 26))
	}
	return 0
}

// mergeChunks sorts chunks by start offset and merges overlapping or adjacent ones.
func mergeChunks(chunks []chunk) []chunk {
	if len(chunks) <= 1 {
		return chunks
	}
	sort.Slice(chunks, func(i, j int) bool { return chunks[i].start < chunks[j].start })
	merged := []chunk{chunks[0]}
	for _, c := range chunks[1:] {
		last := &merged[len(merged)-1]
		if c.start <= last.end {
			if c.end > last.end {
				last.end = c.end
			}
		} else {
			merged = append(merged, c)
		}
	}
	return merged
}

// HasIndex checks if a BAM index file exists.
func HasIndex(bamPath string) bool {
	_, err := os.Stat(bamPath + ".bai")
	return err == nil
}

// BuildIndex creates a BAI index file for a BAM file using real virtual offsets.
func BuildIndex(bamPath string) error {
	f, err := os.Open(bamPath)
	if err != nil {
		return fmt.Errorf("failed to open BAM file: %w", err)
	}
	defer f.Close()

	reader, err := NewReader(f)
	if err != nil {
		return fmt.Errorf("failed to create BAM reader: %w", err)
	}

	header := reader.Header()
	nRef := len(header.References)

	refs := make([]*perRefIndex, nRef)
	for i := range refs {
		refs[i] = &perRefIndex{
			bins:    make(map[uint32][]chunk),
			minVoff: math.MaxInt64,
		}
	}
	var nNoCoor uint64
	var previousRecord *Record

	for {
		virtualOffsetBegin := reader.VirtualOffset()
		record, readErr := reader.Read()
		if readErr != nil {
			if readErr == io.EOF {
				break
			}
			return fmt.Errorf("failed to read BAM record while indexing: %w", readErr)
		}
		virtualOffsetEnd := reader.VirtualOffset()

		if err := validateSortCoordinates(record, nRef); err != nil {
			return fmt.Errorf("cannot index BAM: %w", err)
		}
		if previousRecord != nil && compareCoordinateRecords(previousRecord, record) > 0 {
			return fmt.Errorf("cannot index BAM: records are not coordinate-sorted near %q", record.Name)
		}
		previousRecord = record

		if !recordHasCoordinate(record) {
			nNoCoor++
			continue
		}

		referenceIndex := refs[record.RefID]
		referenceIndex.hasAny = true
		if virtualOffsetBegin < referenceIndex.minVoff {
			referenceIndex.minVoff = virtualOffsetBegin
		}
		if virtualOffsetEnd > referenceIndex.maxVoff {
			referenceIndex.maxVoff = virtualOffsetEnd
		}

		if record.Flags&FlagUnmapped != 0 {
			referenceIndex.nUnmapped++
			continue
		}
		referenceIndex.nMapped++

		begin, end, err := baiReferenceInterval(record)
		if err != nil {
			return fmt.Errorf("cannot index record %q: %w", record.Name, err)
		}

		bin := reg2bin(begin, end)
		referenceIndex.bins[bin] = append(
			referenceIndex.bins[bin],
			chunk{start: virtualOffsetBegin, end: virtualOffsetEnd},
		)

		windowBegin := begin >> 14
		windowEnd := (end - 1) >> 14
		for len(referenceIndex.linear) <= windowEnd {
			referenceIndex.linear = append(referenceIndex.linear, 0)
		}
		for windowIndex := windowBegin; windowIndex <= windowEnd; windowIndex++ {
			if referenceIndex.linear[windowIndex] == 0 || virtualOffsetBegin < referenceIndex.linear[windowIndex] {
				referenceIndex.linear[windowIndex] = virtualOffsetBegin
			}
		}
	}

	// Forward-fill the linear index: empty windows inherit the minimum voff of
	// the preceding non-empty window (optimises region queries).
	for _, ref := range refs {
		for i := 1; i < len(ref.linear); i++ {
			if ref.linear[i] == 0 && ref.linear[i-1] != 0 {
				ref.linear[i] = ref.linear[i-1]
			}
		}
	}

	return writeBAI(bamPath+".bai", nRef, refs, nNoCoor)
}

// writeBAI serialises the index to disk following the BAI format (SAM spec §5.2).
func writeBAI(path string, referenceCount int, references []*perRefIndex, noCoordinateCount uint64) (returnErr error) {
	temporaryFile, err := os.CreateTemp(filepath.Dir(path), filepath.Base(path)+".tmp-*")
	if err != nil {
		return fmt.Errorf("failed to create temporary index file: %w", err)
	}
	temporaryPath := temporaryFile.Name()
	temporaryFileClosed := false
	defer func() {
		if !temporaryFileClosed {
			if closeErr := temporaryFile.Close(); returnErr == nil && closeErr != nil {
				returnErr = fmt.Errorf("failed to close index file: %w", closeErr)
			}
		}
		if returnErr != nil {
			_ = os.Remove(temporaryPath)
		}
	}()

	bufferedWriter := bufio.NewWriterSize(temporaryFile, 1<<20)
	writeUint32 := func(value uint32) error {
		return binary.Write(bufferedWriter, binary.LittleEndian, value)
	}
	writeUint64 := func(value uint64) error {
		return binary.Write(bufferedWriter, binary.LittleEndian, value)
	}

	if _, err := bufferedWriter.Write(baiMagic[:]); err != nil {
		return fmt.Errorf("failed to write BAI magic: %w", err)
	}
	if err := writeUint32(uint32(referenceCount)); err != nil {
		return fmt.Errorf("failed to write BAI reference count: %w", err)
	}

	for _, referenceIndex := range references {
		if !referenceIndex.hasAny {
			if err := writeUint32(0); err != nil {
				return err
			}
			if err := writeUint32(0); err != nil {
				return err
			}
			continue
		}

		if err := writeUint32(uint32(len(referenceIndex.bins)) + 1); err != nil {
			return err
		}
		binIdentifiers := make([]uint32, 0, len(referenceIndex.bins))
		for binIdentifier := range referenceIndex.bins {
			binIdentifiers = append(binIdentifiers, binIdentifier)
		}
		sort.Slice(binIdentifiers, func(leftIndex, rightIndex int) bool {
			return binIdentifiers[leftIndex] < binIdentifiers[rightIndex]
		})
		for _, binIdentifier := range binIdentifiers {
			mergedChunks := mergeChunks(referenceIndex.bins[binIdentifier])
			if err := writeUint32(binIdentifier); err != nil {
				return err
			}
			if err := writeUint32(uint32(len(mergedChunks))); err != nil {
				return err
			}
			for _, indexChunk := range mergedChunks {
				if err := writeUint64(uint64(indexChunk.start)); err != nil {
					return err
				}
				if err := writeUint64(uint64(indexChunk.end)); err != nil {
					return err
				}
			}
		}

		if err := writeUint32(metaBin); err != nil {
			return err
		}
		if err := writeUint32(2); err != nil {
			return err
		}
		minimumVirtualOffset := referenceIndex.minVoff
		if minimumVirtualOffset == math.MaxInt64 {
			minimumVirtualOffset = 0
		}
		for _, value := range []uint64{
			uint64(minimumVirtualOffset),
			uint64(referenceIndex.maxVoff),
			referenceIndex.nMapped,
			referenceIndex.nUnmapped,
		} {
			if err := writeUint64(value); err != nil {
				return err
			}
		}

		if err := writeUint32(uint32(len(referenceIndex.linear))); err != nil {
			return err
		}
		for _, virtualOffset := range referenceIndex.linear {
			if err := writeUint64(uint64(virtualOffset)); err != nil {
				return err
			}
		}
	}

	if err := writeUint64(noCoordinateCount); err != nil {
		return err
	}
	if err := bufferedWriter.Flush(); err != nil {
		return fmt.Errorf("failed to flush index file: %w", err)
	}
	if err := temporaryFile.Sync(); err != nil {
		return fmt.Errorf("failed to sync index file: %w", err)
	}
	if err := temporaryFile.Close(); err != nil {
		return fmt.Errorf("failed to close index file: %w", err)
	}
	temporaryFileClosed = true
	if err := os.Rename(temporaryPath, path); err != nil {
		return fmt.Errorf("failed to install index file: %w", err)
	}
	return nil
}

// EnsureIndex ensures that a BAI index exists for the given BAM file.
func EnsureIndex(bamPath string) error {
	if HasIndex(bamPath) {
		return nil
	}
	return BuildIndex(bamPath)
}
