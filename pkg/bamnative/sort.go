package bamnative

import (
	"fmt"
	"io"
	"os"
	"sort"
)

// SortOptions contains options for BAM sorting
type SortOptions struct {
	OutputPath string // Output path for sorted BAM
	ByName     bool   // Sort by query name instead of coordinate
}

// Sort sorts a BAM file by coordinate
func Sort(inputPath string, opts *SortOptions) error {
	// Open input file
	f, err := os.Open(inputPath)
	if err != nil {
		return fmt.Errorf("failed to open input file: %w", err)
	}
	defer f.Close()

	// Create reader
	reader, err := NewReader(f)
	if err != nil {
		return fmt.Errorf("failed to create reader: %w", err)
	}

	// Read all records
	var records []*Record
	for {
		rec, err := reader.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			return fmt.Errorf("failed to read record: %w", err)
		}
		records = append(records, rec)
	}

	if len(records) == 0 {
		return fmt.Errorf("no records to sort")
	}

	// Sort records by RefID then Position
	sort.Slice(records, func(i, j int) bool {
		if records[i].RefID != records[j].RefID {
			return records[i].RefID < records[j].RefID
		}
		return records[i].Pos < records[j].Pos
	})

	// Create output writer
	header := reader.Header()
	// Update sort order in header
	header.SortOrder = "coordinate"
	if opts.ByName {
		// Sort by query name instead
		sort.Slice(records, func(i, j int) bool {
			return records[i].Name < records[j].Name
		})
		header.SortOrder = "queryname"
	}
	writer, err := NewWriter(opts.OutputPath, header)
	if err != nil {
		return fmt.Errorf("failed to create writer: %w", err)
	}
	defer writer.Close()

	// Write sorted records
	for _, rec := range records {
		if err := writer.Write(rec); err != nil {
			return fmt.Errorf("failed to write record: %w", err)
		}
	}

	return nil
}

// IsSorted checks if a BAM file is coordinate-sorted
func IsSorted(path string) (bool, error) {
	f, err := os.Open(path)
	if err != nil {
		return false, fmt.Errorf("failed to open file: %w", err)
	}
	defer f.Close()

	reader, err := NewReader(f)
	if err != nil {
		return false, fmt.Errorf("failed to create reader: %w", err)
	}

	header := reader.Header()
	return header.SortOrder == "coordinate", nil
}

// sortAndIndexIfNeeded sorts and indexes a BAM file if not already sorted
func SortAndIndexIfNeeded(inputPath, outputPath string) error {
	// Check if already sorted
	isSorted, err := IsSorted(inputPath)
	if err != nil {
		return err
	}

	if isSorted {
		// Just copy the file (or create .bai index)
		return copyBamFile(inputPath, outputPath)
	}

	// Sort the file
	opts := &SortOptions{
		OutputPath: outputPath,
	}
	return Sort(inputPath, opts)
}

// copyBamFile copies a BAM file (used when already sorted)
func copyBamFile(src, dst string) error {
	srcFile, err := os.Open(src)
	if err != nil {
		return err
	}
	defer srcFile.Close()

	dstFile, err := os.Create(dst)
	if err != nil {
		return err
	}
	defer dstFile.Close()

	_, err = io.Copy(dstFile, srcFile)
	return err
}
