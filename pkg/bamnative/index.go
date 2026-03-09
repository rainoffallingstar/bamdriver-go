package bamnative

import (
	"bufio"
	"encoding/binary"
	"fmt"
	"io"
	"math"
	"os"
	"sort"
)

// BAI magic number
var baiMagic = [4]byte{'B', 'A', 'I', 1}

// metaBin is the special BAI pseudo-bin that stores per-reference statistics.
const metaBin = uint32(37450)

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
func alignedLen(rec *Record) int {
	n := 0
	for _, op := range rec.Cigar {
		switch op.Op {
		case CigarMatch, CigarDeletion, CigarSkip, CigarEqual, CigarMismatch:
			n += op.Len
		}
	}
	return n
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

	for {
		voffBeg := reader.VirtualOffset() // virtual offset before reading the record
		rec, err := reader.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			return fmt.Errorf("failed to read BAM record while indexing: %w", err)
		}
		voffEnd := reader.VirtualOffset() // virtual offset after reading the record

		if rec.RefID < 0 || int(rec.RefID) >= nRef {
			nNoCoor++
			continue
		}

		ref := refs[rec.RefID]
		ref.hasAny = true

		// Track min/max voff for the meta bin.
		if voffBeg < ref.minVoff {
			ref.minVoff = voffBeg
		}
		if voffEnd > ref.maxVoff {
			ref.maxVoff = voffEnd
		}

		if rec.Flags&FlagUnmapped != 0 {
			ref.nUnmapped++
			continue // unmapped reads are not added to bin or linear index
		}
		ref.nMapped++

		beg := int(rec.Pos)
		alen := alignedLen(rec)
		end := beg + alen
		if alen == 0 {
			end = beg + 1
		}

		// Bin index
		bin := reg2bin(beg, end)
		ref.bins[bin] = append(ref.bins[bin], chunk{voffBeg, voffEnd})

		// Linear index: one entry per 16 kb window.
		winBeg := beg >> 14
		winEnd := (end - 1) >> 14
		for len(ref.linear) <= winEnd {
			ref.linear = append(ref.linear, 0)
		}
		for w := winBeg; w <= winEnd; w++ {
			if ref.linear[w] == 0 || voffBeg < ref.linear[w] {
				ref.linear[w] = voffBeg
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
func writeBAI(path string, nRef int, refs []*perRefIndex, nNoCoor uint64) error {
	f, err := os.Create(path)
	if err != nil {
		return fmt.Errorf("failed to create index file: %w", err)
	}
	defer f.Close()

	w := bufio.NewWriterSize(f, 1<<20) // 1 MB write buffer

	// Reusable encoding buffers
	var buf4 [4]byte
	var buf8 [8]byte

	write4 := func(v uint32) {
		binary.LittleEndian.PutUint32(buf4[:], v)
		w.Write(buf4[:])
	}
	write8 := func(v uint64) {
		binary.LittleEndian.PutUint64(buf8[:], v)
		w.Write(buf8[:])
	}

	// Magic
	w.Write(baiMagic[:])

	// n_ref
	write4(uint32(nRef))

	for _, ref := range refs {
		if !ref.hasAny {
			// No records for this reference: write empty bin and linear sections.
			write4(0) // n_bin
			write4(0) // n_intv
			continue
		}

		// n_bin = regular bins + 1 meta bin
		write4(uint32(len(ref.bins)) + 1)

		// Write each regular bin (chunks merged to reduce file size).
		for binID, chunks := range ref.bins {
			merged := mergeChunks(chunks)
			write4(binID)
			write4(uint32(len(merged)))
			for _, c := range merged {
				write8(uint64(c.start))
				write8(uint64(c.end))
			}
		}

		// Write meta bin 37450 with 2 pseudo-chunks.
		write4(metaBin)
		write4(2) // n_chunk = 2

		// Pseudo-chunk 1: file range covered by this reference.
		minVoff := ref.minVoff
		if minVoff == math.MaxInt64 {
			minVoff = 0
		}
		write8(uint64(minVoff))
		write8(uint64(ref.maxVoff))

		// Pseudo-chunk 2: mapped / unmapped read counts.
		write8(ref.nMapped)
		write8(ref.nUnmapped)

		// Linear index
		write4(uint32(len(ref.linear)))
		for _, voff := range ref.linear {
			write8(uint64(voff))
		}
	}

	// n_no_coor: reads with no coordinate
	write8(nNoCoor)

	return w.Flush()
}

// EnsureIndex ensures that a BAI index exists for the given BAM file.
func EnsureIndex(bamPath string) error {
	if HasIndex(bamPath) {
		return nil
	}
	return BuildIndex(bamPath)
}
