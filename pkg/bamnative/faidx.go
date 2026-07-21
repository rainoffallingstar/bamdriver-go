package bamnative

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"sync"
)

// FastaIndex represents a FASTA file loaded into memory (backward compatibility)
type FastaIndex struct {
	sequences map[string][]byte
}

// GetSequence returns the sequence for a given reference name
func (f *FastaIndex) GetSequence(refName string) ([]byte, bool) {
	seq, ok := f.sequences[refName]
	return seq, ok
}

// HasSequence checks if a sequence exists
func (f *FastaIndex) HasSequence(refName string) bool {
	_, ok := f.sequences[refName]
	return ok
}

// isGzipped checks if a file is gzip compressed
func isGzipped(path string) bool {
	file, err := os.Open(path)
	if err != nil {
		return false
	}
	defer file.Close()

	// Read first 2 bytes to check for gzip magic number
	buf := make([]byte, 2)
	n, err := file.Read(buf)
	if n < 2 || err != nil {
		return false
	}
	return buf[0] == 0x1f && buf[1] == 0x8b
}

// openFastaFile opens a FASTA file, handling gzip compression
func openFastaFile(path string) (reader io.Reader, closer func() error, err error) {
	if isGzipped(path) {
		file, err := os.Open(path)
		if err != nil {
			return nil, nil, err
		}
		gz, err := gzip.NewReader(file)
		if err != nil {
			file.Close()
			return nil, nil, err
		}
		return gz, func() error {
			gz.Close()
			return file.Close()
		}, nil
	}
	f, err := os.Open(path)
	return f, f.Close, err
}

// gzipReader wraps gzip reader with proper cleanup
type gzipReader struct {
	file   *os.File
	reader *gzip.Reader
}

func (r *gzipReader) Read(p []byte) (n int, err error) {
	return r.reader.Read(p)
}

func (r *gzipReader) Close() error {
	r.reader.Close()
	return r.file.Close()
}

// FastaIndexEntry represents an entry in the FASTA index
type FastaIndexEntry struct {
	Name   string
	Length int64
	Offset int64
	LineB  int64
	LineL  int64
}

// FastaReader provides streaming FASTA access with optional indexing
type FastaReader struct {
	path       string
	isGzipped  bool
	entries    map[string]*FastaIndexEntry
	entryOrder []string
	cache      map[string][]byte
	maxCache   int
	mutex      sync.RWMutex
}

// NewFastaReader creates a new FASTA reader with optional indexing
func NewFastaReader(path string) (*FastaReader, error) {
	gzipped := isGzipped(path)

	reader := &FastaReader{
		path:      path,
		isGzipped: gzipped,
		entries:   make(map[string]*FastaIndexEntry),
		cache:     make(map[string][]byte),
		maxCache:  10,
	}

	if gzipped {
		entries, err := reader.buildIndexGzipped()
		if err != nil {
			return nil, err
		}
		reader.entries = entries
		return reader, nil
	}

	indexPath := path + ".fai"
	if _, err := os.Stat(indexPath); err == nil {
		if err := reader.loadIndex(indexPath); err != nil {
			if rebuildErr := reader.buildAndSaveIndex(indexPath); rebuildErr != nil {
				return nil, fmt.Errorf("invalid FASTA index (%v), and rebuilding failed: %w", err, rebuildErr)
			}
		}
	} else if os.IsNotExist(err) {
		if err := reader.buildAndSaveIndex(indexPath); err != nil {
			return nil, err
		}
	} else {
		return nil, fmt.Errorf("failed to inspect FASTA index: %w", err)
	}

	return reader, nil
}

// buildAndSaveIndex builds an index and installs it atomically.
func (fr *FastaReader) buildAndSaveIndex(indexPath string) (returnErr error) {
	entries, err := fr.buildIndex()
	if err != nil {
		return err
	}

	temporaryFile, err := os.CreateTemp(filepath.Dir(indexPath), filepath.Base(indexPath)+".tmp-*")
	if err != nil {
		return fmt.Errorf("failed to create temporary FASTA index: %w", err)
	}
	temporaryPath := temporaryFile.Name()
	temporaryFileClosed := false
	defer func() {
		if !temporaryFileClosed {
			if closeErr := temporaryFile.Close(); returnErr == nil && closeErr != nil {
				returnErr = closeErr
			}
		}
		if returnErr != nil {
			_ = os.Remove(temporaryPath)
		}
	}()

	bufferedWriter := bufio.NewWriter(temporaryFile)
	for _, name := range fr.entryOrder {
		entry := entries[name]
		if _, err := fmt.Fprintf(
			bufferedWriter,
			"%s\t%d\t%d\t%d\t%d\n",
			entry.Name,
			entry.Length,
			entry.Offset,
			entry.LineB,
			entry.LineL,
		); err != nil {
			return fmt.Errorf("failed to write FASTA index: %w", err)
		}
	}
	if err := bufferedWriter.Flush(); err != nil {
		return fmt.Errorf("failed to flush FASTA index: %w", err)
	}
	if err := temporaryFile.Sync(); err != nil {
		return fmt.Errorf("failed to sync FASTA index: %w", err)
	}
	if err := temporaryFile.Close(); err != nil {
		return fmt.Errorf("failed to close FASTA index: %w", err)
	}
	temporaryFileClosed = true
	if err := os.Rename(temporaryPath, indexPath); err != nil {
		return fmt.Errorf("failed to install FASTA index: %w", err)
	}
	fr.entries = entries
	return nil
}

// buildIndex scans physical FASTA lines so offsets include their real line endings.
func (fr *FastaReader) buildIndex() (map[string]*FastaIndexEntry, error) {
	if fr.isGzipped {
		return fr.buildIndexGzipped()
	}

	file, err := os.Open(fr.path)
	if err != nil {
		return nil, fmt.Errorf("failed to open FASTA file: %w", err)
	}
	defer file.Close()

	entries := make(map[string]*FastaIndexEntry)
	fr.entryOrder = nil
	bufferedReader := bufio.NewReader(file)
	var currentEntry *FastaIndexEntry
	var offset int64
	var sawShortSequenceLine bool

	for {
		physicalLine, readErr := bufferedReader.ReadBytes('\n')
		if len(physicalLine) > 0 {
			lineStart := offset
			offset += int64(len(physicalLine))
			sequenceLine := bytes.TrimSuffix(physicalLine, []byte{'\n'})
			sequenceLine = bytes.TrimSuffix(sequenceLine, []byte{'\r'})

			if len(sequenceLine) == 0 {
				return nil, fmt.Errorf("blank lines are not supported in indexed FASTA")
			}
			if sequenceLine[0] == '>' {
				name := parseFastaName(string(sequenceLine[1:]))
				if name == "" {
					return nil, fmt.Errorf("FASTA header at offset %d has no sequence name", lineStart)
				}
				if _, exists := entries[name]; exists {
					return nil, fmt.Errorf("duplicate FASTA sequence name %q", name)
				}
				currentEntry = &FastaIndexEntry{Name: name, Offset: offset}
				entries[name] = currentEntry
				fr.entryOrder = append(fr.entryOrder, name)
				sawShortSequenceLine = false
			} else {
				if currentEntry == nil {
					return nil, fmt.Errorf("sequence data appears before the first FASTA header")
				}
				lineBases := int64(len(sequenceLine))
				lineWidth := int64(len(physicalLine))
				lineHasTerminator := physicalLine[len(physicalLine)-1] == '\n'
				if lineBases == 0 {
					return nil, fmt.Errorf("empty sequence line for %q", currentEntry.Name)
				}
				if currentEntry.LineB == 0 {
					currentEntry.LineB = lineBases
					currentEntry.LineL = lineWidth
				} else {
					if sawShortSequenceLine {
						return nil, fmt.Errorf("sequence %q has data after its final short line", currentEntry.Name)
					}
					if lineBases > currentEntry.LineB || (lineBases == currentEntry.LineB && lineWidth != currentEntry.LineL && lineHasTerminator) {
						return nil, fmt.Errorf("sequence %q has inconsistent line widths", currentEntry.Name)
					}
					if lineBases < currentEntry.LineB || !lineHasTerminator {
						sawShortSequenceLine = true
					}
				}
				currentEntry.Length += lineBases
			}
		}
		if readErr != nil {
			if readErr == io.EOF {
				break
			}
			return nil, fmt.Errorf("failed to scan FASTA: %w", readErr)
		}
	}

	if len(entries) == 0 {
		return nil, fmt.Errorf("FASTA contains no sequences")
	}
	return entries, nil
}

// buildIndexGzipped builds index for gzipped FASTA by loading into memory.
// For gzipped files we cannot seek, so we load all sequences into memory in a
// single pass and store them directly into fr.cache.  Previous versions ran N
// additional full-file scans (one per chromosome) after building the index,
// making initialisation O(N²) – fixed here.
func (fr *FastaReader) buildIndexGzipped() (map[string]*FastaIndexEntry, error) {
	reader, closer, err := openFastaFile(fr.path)
	if err != nil {
		return nil, fmt.Errorf("failed to open FASTA file: %w", err)
	}
	defer closer()

	entries := make(map[string]*FastaIndexEntry)
	var currentName string
	var currentSeq []byte

	scanner := bufio.NewScanner(reader)
	// Increase buffer size for long lines
	buf := make([]byte, 0, 1024*1024)
	scanner.Buffer(buf, 1024*1024)

	for scanner.Scan() {
		line := scanner.Text()
		line = strings.TrimSpace(line)

		if len(line) == 0 {
			continue
		}

		if line[0] == '>' {
			// Flush previous sequence directly into cache (no second scan needed)
			if currentName != "" {
				entries[currentName] = &FastaIndexEntry{
					Name:   currentName,
					Length: int64(len(currentSeq)),
					LineL:  0, // Not used for gzipped
				}
				fr.cache[currentName] = currentSeq
			}

			// Parse new sequence name
			currentName = strings.TrimSpace(line[1:])
			if idx := strings.Index(currentName, " "); idx != -1 {
				currentName = currentName[:idx]
			}
			currentSeq = nil
		} else {
			currentSeq = append(currentSeq, line...)
		}
	}

	// Flush last sequence
	if currentName != "" {
		entries[currentName] = &FastaIndexEntry{
			Name:   currentName,
			Length: int64(len(currentSeq)),
			LineL:  0,
		}
		fr.cache[currentName] = currentSeq
	}

	if err := scanner.Err(); err != nil {
		return nil, err
	}

	return entries, nil
}

// loadSequenceGzipped loads a sequence from gzipped FASTA
func (fr *FastaReader) loadSequenceGzipped(name string) ([]byte, bool) {
	// Check cache first
	if seq, ok := fr.cache[name]; ok {
		return seq, true
	}

	// For gzipped files, we need to scan from beginning
	reader, closer, err := openFastaFile(fr.path)
	if err != nil {
		return nil, false
	}
	defer closer()

	var currentName string
	var currentSeq []byte

	scanner := bufio.NewScanner(reader)
	buf := make([]byte, 0, 1024*1024)
	scanner.Buffer(buf, 1024*1024)

	for scanner.Scan() {
		line := scanner.Text()
		line = strings.TrimSpace(line)

		if len(line) == 0 {
			continue
		}

		if line[0] == '>' {
			if currentName == name {
				return currentSeq, len(currentSeq) > 0
			}
			currentName = strings.TrimSpace(line[1:])
			if idx := strings.Index(currentName, " "); idx != -1 {
				currentName = currentName[:idx]
			}
			currentSeq = nil
		} else {
			currentSeq = append(currentSeq, line...)
		}
	}

	// Check last sequence
	if currentName == name {
		return currentSeq, len(currentSeq) > 0
	}

	return nil, false
}

// loadIndex loads and validates a FASTA index file.
func (fr *FastaReader) loadIndex(indexPath string) error {
	file, err := os.Open(indexPath)
	if err != nil {
		return fmt.Errorf("failed to open index file: %w", err)
	}
	defer file.Close()

	entries := make(map[string]*FastaIndexEntry)
	entryOrder := make([]string, 0)
	scanner := bufio.NewScanner(file)
	lineNumber := 0
	for scanner.Scan() {
		lineNumber++
		fields := strings.Split(scanner.Text(), "\t")
		if len(fields) != 5 || fields[0] == "" {
			return fmt.Errorf("invalid FASTA index line %d", lineNumber)
		}
		if _, exists := entries[fields[0]]; exists {
			return fmt.Errorf("duplicate FASTA index entry %q", fields[0])
		}
		values := make([]int64, 4)
		for valueIndex := range values {
			parsedValue, parseErr := strconv.ParseInt(fields[valueIndex+1], 10, 64)
			if parseErr != nil || parsedValue < 0 {
				return fmt.Errorf("invalid FASTA index value on line %d", lineNumber)
			}
			values[valueIndex] = parsedValue
		}
		entry := &FastaIndexEntry{
			Name:   fields[0],
			Length: values[0],
			Offset: values[1],
			LineB:  values[2],
			LineL:  values[3],
		}
		if entry.Length > 0 && (entry.LineB <= 0 || entry.LineL < entry.LineB) {
			return fmt.Errorf("invalid FASTA line geometry for %q", entry.Name)
		}
		entries[entry.Name] = entry
		entryOrder = append(entryOrder, entry.Name)
	}
	if err := scanner.Err(); err != nil {
		return fmt.Errorf("failed to read FASTA index: %w", err)
	}
	if len(entries) == 0 {
		return fmt.Errorf("FASTA index is empty")
	}
	fr.entries = entries
	fr.entryOrder = entryOrder
	return nil
}

func parseFastaName(headerText string) string {
	fields := strings.Fields(headerText)
	if len(fields) == 0 {
		return ""
	}
	return fields[0]
}

// HasIndex returns true if index is available.
func (fr *FastaReader) HasIndex() bool {
	fr.mutex.RLock()
	defer fr.mutex.RUnlock()
	return len(fr.entries) > 0
}

// GetSequence retrieves a sequence by name using a concurrency-safe cache.
func (fr *FastaReader) GetSequence(name string) ([]byte, bool) {
	fr.mutex.Lock()
	defer fr.mutex.Unlock()

	if sequence, ok := fr.cache[name]; ok {
		return sequence, true
	}

	// Try to load from file
	var seq []byte
	var ok bool

	// For gzipped files, always use gzipped loader
	if fr.isGzipped {
		seq, ok = fr.loadSequenceGzipped(name)
		if !ok {
			// Two-directional "chr" prefix fallback so that mismatches between
			// BAM header chromosome names and FASTA names are resolved:
			//   "chr1" (BAM) + FASTA has "1"  → strip prefix
			//   "1"    (BAM) + FASTA has "chr1" → add prefix
			if strings.HasPrefix(name, "chr") {
				altName := name[3:] // "chr1" → "1"
				if seq, ok = fr.loadSequenceGzipped(altName); !ok {
					altName = "chr" + name // last resort: "chrchr1"
					seq, ok = fr.loadSequenceGzipped(altName)
				}
			} else {
				altName := "chr" + name // "1" → "chr1"
				seq, ok = fr.loadSequenceGzipped(altName)
			}
		}
	} else if fr.entries != nil {
		// Use index to seek to sequence
		seq, ok = fr.loadSequenceByIndex(name)
	} else {
		// Scan from beginning (slower)
		seq, ok = fr.loadSequenceByScan(name)
	}

	if !ok {
		return nil, false
	}

	// Add to cache with LRU eviction
	fr.addToCache(name, seq)

	return seq, true
}

// loadSequenceByIndex loads sequence using index file
func (fr *FastaReader) loadSequenceByIndex(name string) ([]byte, bool) {
	entry, ok := fr.entries[name]
	if !ok {
		// Two-directional "chr" prefix fallback:
		//   "chr1" (BAM) + FASTA has "1"   → strip prefix
		//   "1"    (BAM) + FASTA has "chr1" → add prefix
		if strings.HasPrefix(name, "chr") {
			altName := name[3:] // "chr1" → "1"
			if entry, ok = fr.entries[altName]; !ok {
				altName = "chr" + name // last resort: "chrchr1"
				entry, ok = fr.entries[altName]
			}
		} else {
			altName := "chr" + name // "1" → "chr1"
			entry, ok = fr.entries[altName]
		}
		if !ok {
			return nil, false
		}
	}

	file, err := os.Open(fr.path)
	if err != nil {
		return nil, false
	}
	defer file.Close()

	if entry.Length == 0 {
		return []byte{}, true
	}
	if entry.LineB <= 0 || entry.LineL < entry.LineB {
		return nil, false
	}
	sequence := make([]byte, 0, entry.Length)
	remainingBases := entry.Length
	fileOffset := entry.Offset
	for remainingBases > 0 {
		basesToRead := entry.LineB
		if remainingBases < basesToRead {
			basesToRead = remainingBases
		}
		lineData := make([]byte, basesToRead)
		if _, err := file.ReadAt(lineData, fileOffset); err != nil {
			return nil, false
		}
		sequence = append(sequence, lineData...)
		remainingBases -= basesToRead
		fileOffset += entry.LineL
	}
	return sequence, true
}

// loadSequenceByScan loads sequence by scanning from beginning
func (fr *FastaReader) loadSequenceByScan(name string) ([]byte, bool) {
	file, err := os.Open(fr.path)
	if err != nil {
		return nil, false
	}
	defer file.Close()

	var currentName string
	var currentSeq []byte

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := scanner.Text()
		line = strings.TrimSpace(line)

		if len(line) == 0 {
			continue
		}

		if line[0] == '>' {
			if currentName == name {
				return currentSeq, len(currentSeq) > 0
			}
			// Parse new name
			currentName = strings.TrimSpace(line[1:])
			if idx := strings.Index(currentName, " "); idx != -1 {
				currentName = currentName[:idx]
			}
			currentSeq = nil
		} else {
			currentSeq = append(currentSeq, line...)
		}
	}

	// Check last sequence
	if currentName == name {
		return currentSeq, len(currentSeq) > 0
	}

	return nil, false
}

// addToCache adds sequence to cache with LRU eviction
func (fr *FastaReader) addToCache(name string, sequence []byte) {
	if _, exists := fr.cache[name]; exists {
		fr.cache[name] = sequence
		return
	}
	if len(fr.cache) >= fr.maxCache {
		for cachedName := range fr.cache {
			delete(fr.cache, cachedName)
			break
		}
	}
	fr.cache[name] = sequence
}

// LoadFasta loads entire FASTA into memory (backward compatibility)
func LoadFasta(path string) (*FastaIndex, error) {
	file, err := os.Open(path)
	if err != nil {
		return nil, fmt.Errorf("failed to open FASTA file: %w", err)
	}
	defer file.Close()

	sequences := make(map[string][]byte)
	var currentName string
	var currentSeq []byte

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := scanner.Text()
		line = strings.TrimSpace(line)

		if len(line) == 0 {
			continue
		}

		if line[0] == '>' {
			if currentName != "" {
				sequences[currentName] = currentSeq
			}
			currentName = strings.TrimSpace(line[1:])
			if idx := strings.Index(currentName, " "); idx != -1 {
				currentName = currentName[:idx]
			}
			currentSeq = nil
		} else {
			currentSeq = append(currentSeq, line...)
		}
	}

	if currentName != "" {
		sequences[currentName] = currentSeq
	}

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("error reading FASTA file: %w", err)
	}

	return &FastaIndex{sequences: sequences}, nil
}
