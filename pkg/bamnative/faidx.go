package bamnative

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"strings"
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
	path      string
	isGzipped bool
	entries   map[string]*FastaIndexEntry
	cache     map[string][]byte
	cacheSize int
	maxCache  int // Maximum number of sequences to keep in cache
}

// NewFastaReader creates a new FASTA reader with optional indexing
func NewFastaReader(path string) (*FastaReader, error) {
	gzipped := isGzipped(path)

	fr := &FastaReader{
		path:      path,
		isGzipped: gzipped,
		entries:   make(map[string]*FastaIndexEntry),
		cache:     make(map[string][]byte),
		cacheSize: 0,
		maxCache:  10, // Keep up to 10 sequences in cache
	}

	// For gzipped files, we cannot use index - load all into memory
	if gzipped {
		fmt.Fprintf(os.Stderr, "Loading gzipped FASTA: %s\n", path)
		entries, err := fr.buildIndexGzipped()
		if err != nil {
			fmt.Fprintf(os.Stderr, "Warning: failed to load gzipped FASTA: %v\n", err)
		} else {
			fr.entries = entries
		}
		return fr, nil
	}

	// Try to load index file
	indexPath := path + ".fai"
	if _, err := os.Stat(indexPath); err == nil {
		if err := fr.loadIndex(indexPath); err != nil {
			// Index loading failed, continue without index
			fr.entries = nil
		}
	} else {
		// Index file not found, generate it
		fmt.Fprintf(os.Stderr, "Generating FASTA index for: %s\n", path)
		if err := fr.buildAndSaveIndex(indexPath); err != nil {
			fmt.Fprintf(os.Stderr, "Warning: failed to generate index: %v\n", err)
			fr.entries = nil
		}
	}

	return fr, nil
}

// buildAndSaveIndex builds index from FASTA and saves to file
func (fr *FastaReader) buildAndSaveIndex(indexPath string) error {
	entries, err := fr.buildIndex()
	if err != nil {
		return err
	}
	fr.entries = entries

	// Save index to file
	file, err := os.Create(indexPath)
	if err != nil {
		return fmt.Errorf("failed to create index file: %w", err)
	}
	defer file.Close()

	for _, entry := range entries {
		fmt.Fprintf(file, "%s\t%d\t%d\t%d\t%d\n",
			entry.Name, entry.Length, entry.Offset, entry.LineB, entry.LineL)
	}

	fmt.Fprintf(os.Stderr, "FASTA index saved to: %s\n", indexPath)
	return nil
}

// buildIndex builds FASTA index by scanning the file
func (fr *FastaReader) buildIndex() (map[string]*FastaIndexEntry, error) {
	// For gzipped files, we cannot build index - load all into memory instead
	if fr.isGzipped {
		return fr.buildIndexGzipped()
	}

	file, err := os.Open(fr.path)
	if err != nil {
		return nil, fmt.Errorf("failed to open FASTA file: %w", err)
	}
	defer file.Close()

	entries := make(map[string]*FastaIndexEntry)

	var currentName string
	var seqStartOffset int64
	var lineB int64 = 0
	var lineL int64 = 0

	scanner := bufio.NewScanner(file)
	offset := int64(0)
	firstSeq := true

	for scanner.Scan() {
		line := scanner.Text()
		lineLen := int64(len(line))

		// Detect line format
		if firstSeq && lineLen > 0 {
			if line[len(line)-1] == '\r' {
				lineL = lineLen - 1 // Windows line ending
			} else {
				lineL = lineLen // Unix line ending
			}
			// Check if lines have fixed length (binary FASTA)
			if line[0] != '>' {
				// This is sequence line, determine lineB
				lineB = lineL
			}
		}

		if len(line) == 0 {
			offset += int64(lineLen) + 1 // +1 for newline
			continue
		}

		if line[0] == '>' {
			// Save previous sequence info
			if !firstSeq && currentName != "" {
				if entry, ok := entries[currentName]; ok {
					entry.Length = offset - seqStartOffset - int64(lineB) - 1
				}
			}

			// Parse new sequence name
			currentName = strings.TrimSpace(line[1:])
			if idx := strings.Index(currentName, " "); idx != -1 {
				currentName = currentName[:idx]
			}

			// Record start of this sequence
			seqStartOffset = offset + int64(len(line)) + 1
			entries[currentName] = &FastaIndexEntry{
				Name:   currentName,
				Offset: seqStartOffset,
				LineB:  lineB,
				LineL:  lineL,
			}
			firstSeq = false
		}

		offset += int64(lineLen) + 1 // +1 for newline
	}

	// Save last sequence length
	if currentName != "" {
		if entry, ok := entries[currentName]; ok {
			// Calculate total file size
			fileInfo, _ := file.Stat()
			entry.Length = fileInfo.Size() - entry.Offset - int64(entry.LineB)
			// Adjust for number of newlines
			if entry.Length > 0 && entry.LineB > 0 {
				entry.Length = entry.Length - (entry.Length / int64(entry.LineB))
			}
		}
	}

	return entries, scanner.Err()
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

// loadIndex loads the FASTA index file
func (fr *FastaReader) loadIndex(indexPath string) error {
	file, err := os.Open(indexPath)
	if err != nil {
		return fmt.Errorf("failed to open index file: %w", err)
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := scanner.Text()
		fields := strings.Fields(line)
		if len(fields) < 5 {
			continue
		}

		entry := &FastaIndexEntry{
			Name:   fields[0],
			Length: 0,
			Offset: 0,
			LineB:  0,
			LineL:  0,
		}
		fmt.Sscanf(fields[1], "%d", &entry.Length)
		fmt.Sscanf(fields[2], "%d", &entry.Offset)
		fmt.Sscanf(fields[3], "%d", &entry.LineB)
		fmt.Sscanf(fields[4], "%d", &entry.LineL)

		fr.entries[entry.Name] = entry
	}

	return scanner.Err()
}

// HasIndex returns true if index is available
func (fr *FastaReader) HasIndex() bool {
	return fr.entries != nil && len(fr.entries) > 0
}

// GetSequence retrieves a sequence by name, using cache
func (fr *FastaReader) GetSequence(name string) ([]byte, bool) {
	// Check cache first
	if seq, ok := fr.cache[name]; ok {
		return seq, true
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

	// Seek to sequence position
	_, err = file.Seek(entry.Offset, 0)
	if err != nil {
		return nil, false
	}

	// Read sequence lines
	var seq []byte
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := scanner.Text()
		if len(line) == 0 {
			continue
		}
		if line[0] == '>' {
			break
		}
		seq = append(seq, line...)
	}

	return seq, len(seq) > 0
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
func (fr *FastaReader) addToCache(name string, seq []byte) {
	// Evict if cache is full
	if len(fr.cache) >= fr.maxCache {
		// Simple eviction: remove first entry
		for k := range fr.cache {
			delete(fr.cache, k)
			break
		}
	}

	fr.cache[name] = seq
	fr.cacheSize++
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
