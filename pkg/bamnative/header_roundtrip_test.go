package bamnative

import (
	"encoding/binary"
	"os"
	"path/filepath"
	"reflect"
	"testing"

	"github.com/rainoffallingstar/bamdriver-go/pkg/bgzip"
)

func TestHeaderRoundTripPreservesVersionAndOrderedDuplicateLines(t *testing.T) {
	temporaryDirectory := t.TempDir()
	bamPath := filepath.Join(temporaryDirectory, "header.bam")
	header := &Header{
		Version:   "1.4",
		SortOrder: "coordinate",
		References: []*Reference{
			{ID: 0, Name: "chr1", Len: 1000},
		},
		OtherHeaderLines: []string{
			"@CO\tfirst comment",
			"@RG\tID:group-one\tSM:sample",
			"@CO\tsecond comment",
			"@PG\tID:program-one\tPN:test",
		},
	}

	writer, err := NewWriter(bamPath, header)
	if err != nil {
		t.Fatalf("NewWriter: %v", err)
	}
	if err := writer.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	file, err := os.Open(bamPath)
	if err != nil {
		t.Fatalf("Open: %v", err)
	}
	defer file.Close()
	reader, err := NewReader(file)
	if err != nil {
		t.Fatalf("NewReader: %v", err)
	}
	if reader.Header().Version != header.Version {
		t.Fatalf("version = %q, want %q", reader.Header().Version, header.Version)
	}
	if !reflect.DeepEqual(reader.Header().OtherHeaderLines, header.OtherHeaderLines) {
		t.Fatalf("header lines = %#v, want %#v", reader.Header().OtherHeaderLines, header.OtherHeaderLines)
	}
}

func TestReaderAcceptsBinaryReferenceDictionaryWithoutTextSQ(t *testing.T) {
	temporaryDirectory := t.TempDir()
	bamPath := filepath.Join(temporaryDirectory, "binary-only-reference.bam")
	writeHeaderOnlyBAM(t, bamPath, "@HD\tVN:1.6\tSO:unknown\n", "chr1", 1000)

	file, err := os.Open(bamPath)
	if err != nil {
		t.Fatalf("Open: %v", err)
	}
	defer file.Close()
	reader, err := NewReader(file)
	if err != nil {
		t.Fatalf("NewReader: %v", err)
	}
	if len(reader.Header().References) != 1 {
		t.Fatalf("reference count = %d, want 1", len(reader.Header().References))
	}
	reference := reader.Header().References[0]
	if reference.Name != "chr1" || reference.Len != 1000 {
		t.Fatalf("reference = %#v", reference)
	}
}

func TestReaderRejectsTextAndBinaryReferenceMismatch(t *testing.T) {
	testCases := []struct {
		name         string
		binaryName   string
		binaryLength int32
	}{
		{name: "name", binaryName: "chr2", binaryLength: 1000},
		{name: "length", binaryName: "chr1", binaryLength: 999},
	}

	for _, testCase := range testCases {
		t.Run(testCase.name, func(t *testing.T) {
			temporaryDirectory := t.TempDir()
			bamPath := filepath.Join(temporaryDirectory, "mismatch.bam")
			writeHeaderOnlyBAM(t, bamPath, "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n", testCase.binaryName, testCase.binaryLength)

			file, err := os.Open(bamPath)
			if err != nil {
				t.Fatalf("Open: %v", err)
			}
			defer file.Close()
			if _, err := NewReader(file); err == nil {
				t.Fatal("NewReader accepted mismatched text and binary reference dictionaries")
			}
		})
	}
}

func writeHeaderOnlyBAM(t *testing.T, path, headerText, referenceName string, referenceLength int32) {
	t.Helper()
	writer, err := bgzip.NewWriter(path)
	if err != nil {
		t.Fatalf("NewWriter: %v", err)
	}
	write := func(value interface{}) {
		t.Helper()
		if err := binary.Write(writer, binary.LittleEndian, value); err != nil {
			_ = writer.Close()
			t.Fatalf("binary.Write: %v", err)
		}
	}
	if _, err := writer.Write([]byte{'B', 'A', 'M', 1}); err != nil {
		t.Fatalf("write magic: %v", err)
	}
	write(int32(len(headerText)))
	if _, err := writer.Write([]byte(headerText)); err != nil {
		t.Fatalf("write header text: %v", err)
	}
	write(int32(1))
	write(int32(len(referenceName) + 1))
	if _, err := writer.Write(append([]byte(referenceName), 0)); err != nil {
		t.Fatalf("write reference name: %v", err)
	}
	write(referenceLength)
	if err := writer.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
}
