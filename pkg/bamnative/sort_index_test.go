package bamnative

import (
	"fmt"
	"io"
	"os"
	"path/filepath"
	"reflect"
	"strings"
	"testing"
)

func TestSortExternalMergePreservesStableCoordinateOrder(t *testing.T) {
	temporaryDirectory := t.TempDir()
	inputPath := filepath.Join(temporaryDirectory, "many-runs-input.bam")
	outputPath := filepath.Join(temporaryDirectory, "many-runs-output.bam")
	header := testCoordinateHeader()
	header.SortOrder = "unknown"

	records := make([]*Record, 0, 70)
	wantNames := make([]string, 0, 70)
	for recordIndex := 0; recordIndex < 70; recordIndex++ {
		name := fmt.Sprintf("record-%03d", recordIndex)
		records = append(records, testMappedRecord(name, 0, 10))
		wantNames = append(wantNames, name)
	}
	writeTestBAM(t, inputPath, header, records)

	if err := Sort(inputPath, &SortOptions{
		OutputPath:         outputPath,
		MemoryLimitBytes:   1,
		TemporaryDirectory: temporaryDirectory,
	}); err != nil {
		t.Fatalf("Sort external merge: %v", err)
	}

	gotNames := readRecordNames(t, outputPath)
	if !reflect.DeepEqual(gotNames, wantNames) {
		t.Fatalf("stable order mismatch:\n got: %#v\nwant: %#v", gotNames, wantNames)
	}
	temporaryMatches, err := filepath.Glob(filepath.Join(temporaryDirectory, ".bam-sort-*"))
	if err != nil {
		t.Fatalf("Glob: %v", err)
	}
	if len(temporaryMatches) != 0 {
		t.Fatalf("temporary sort directories were not removed: %v", temporaryMatches)
	}
}

func TestIsSortedScansRecordsInsteadOfTrustingHeader(t *testing.T) {
	temporaryDirectory := t.TempDir()
	bamPath := filepath.Join(temporaryDirectory, "declared-sorted.bam")
	header := testCoordinateHeader()
	writeTestBAM(t, bamPath, header, []*Record{
		testMappedRecord("later", 0, 20),
		testMappedRecord("earlier", 0, 10),
	})

	isSorted, err := IsSorted(bamPath)
	if err != nil {
		t.Fatalf("IsSorted: %v", err)
	}
	if isSorted {
		t.Fatal("IsSorted returned true for records in descending coordinate order")
	}
}

func TestSortPlacesRecordsWithoutCoordinatesLast(t *testing.T) {
	temporaryDirectory := t.TempDir()
	inputPath := filepath.Join(temporaryDirectory, "input.bam")
	outputPath := filepath.Join(temporaryDirectory, "sorted.bam")
	header := testCoordinateHeader()
	header.SortOrder = "unknown"
	writeTestBAM(t, inputPath, header, []*Record{
		testUnplacedRecord("unplaced"),
		testMappedRecord("second-reference", 1, 1),
		testMappedRecord("later", 0, 20),
		testMappedRecord("earlier", 0, 10),
	})

	if err := Sort(inputPath, &SortOptions{OutputPath: outputPath}); err != nil {
		t.Fatalf("Sort: %v", err)
	}

	names := readRecordNames(t, outputPath)
	wantNames := []string{"earlier", "later", "second-reference", "unplaced"}
	if len(names) != len(wantNames) {
		t.Fatalf("record count = %d, want %d", len(names), len(wantNames))
	}
	for recordIndex := range wantNames {
		if names[recordIndex] != wantNames[recordIndex] {
			t.Fatalf("record %d = %q, want %q", recordIndex, names[recordIndex], wantNames[recordIndex])
		}
	}

	isSorted, err := IsSorted(outputPath)
	if err != nil {
		t.Fatalf("IsSorted sorted output: %v", err)
	}
	if !isSorted {
		t.Fatal("sorted output was not recognized as coordinate-sorted")
	}
}

func TestBuildIndexRejectsUnsortedRecords(t *testing.T) {
	temporaryDirectory := t.TempDir()
	bamPath := filepath.Join(temporaryDirectory, "unsorted.bam")
	writeTestBAM(t, bamPath, testCoordinateHeader(), []*Record{
		testMappedRecord("later", 0, 20),
		testMappedRecord("earlier", 0, 10),
	})

	err := BuildIndex(bamPath)
	if err == nil {
		t.Fatal("BuildIndex succeeded for unsorted records")
	}
	if !strings.Contains(err.Error(), "not coordinate-sorted") {
		t.Fatalf("BuildIndex error = %v, want coordinate-sort error", err)
	}
	if _, statErr := os.Stat(bamPath + ".bai"); !os.IsNotExist(statErr) {
		t.Fatalf("unexpected index file after failed BuildIndex: %v", statErr)
	}
}

func TestSortAndIndexIfNeededCreatesSortedBAMAndIndex(t *testing.T) {
	temporaryDirectory := t.TempDir()
	inputPath := filepath.Join(temporaryDirectory, "input.bam")
	outputPath := filepath.Join(temporaryDirectory, "output.bam")
	header := testCoordinateHeader()
	header.SortOrder = "unknown"
	writeTestBAM(t, inputPath, header, []*Record{
		testUnplacedRecord("unplaced"),
		testMappedRecord("later", 0, 20),
		testMappedRecord("earlier", 0, 10),
	})

	if err := SortAndIndexIfNeeded(inputPath, outputPath); err != nil {
		t.Fatalf("SortAndIndexIfNeeded: %v", err)
	}
	if !HasIndex(outputPath) {
		t.Fatal("SortAndIndexIfNeeded did not create a BAI index")
	}
	indexData, err := os.ReadFile(outputPath + ".bai")
	if err != nil {
		t.Fatalf("ReadFile index: %v", err)
	}
	if len(indexData) < len(baiMagic) || string(indexData[:len(baiMagic)]) != string(baiMagic[:]) {
		t.Fatalf("index does not begin with BAI magic: %v", indexData)
	}
	isSorted, err := IsSorted(outputPath)
	if err != nil {
		t.Fatalf("IsSorted output: %v", err)
	}
	if !isSorted {
		t.Fatal("SortAndIndexIfNeeded output is not sorted")
	}
}

func testCoordinateHeader() *Header {
	return &Header{
		Version:   "1.6",
		SortOrder: "coordinate",
		References: []*Reference{
			{ID: 0, Name: "chr1", Len: 1000},
			{ID: 1, Name: "chr2", Len: 1000},
		},
	}
}

func testMappedRecord(name string, referenceID, position int32) *Record {
	return &Record{
		Name:      name,
		RefID:     referenceID,
		Pos:       position,
		MateRefID: -1,
		MatePos:   -1,
		Cigar:     []CigarOp{{Op: CigarMatch, Len: 1}},
		Seq:       "A",
	}
}

func testUnplacedRecord(name string) *Record {
	return &Record{
		Name:      name,
		Flags:     FlagUnmapped,
		RefID:     -1,
		Pos:       -1,
		MateRefID: -1,
		MatePos:   -1,
	}
}

func writeTestBAM(t *testing.T, path string, header *Header, records []*Record) {
	t.Helper()
	writer, err := NewWriter(path, header)
	if err != nil {
		t.Fatalf("NewWriter: %v", err)
	}
	for _, record := range records {
		if err := writer.Write(record); err != nil {
			_ = writer.Close()
			t.Fatalf("Write %q: %v", record.Name, err)
		}
	}
	if err := writer.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
}

func readRecordNames(t *testing.T, path string) []string {
	t.Helper()
	file, err := os.Open(path)
	if err != nil {
		t.Fatalf("Open: %v", err)
	}
	defer file.Close()
	reader, err := NewReader(file)
	if err != nil {
		t.Fatalf("NewReader: %v", err)
	}

	var names []string
	for {
		record, readErr := reader.Read()
		if readErr != nil {
			if readErr == io.EOF {
				return names
			}
			t.Fatalf("Read: %v", readErr)
		}
		names = append(names, record.Name)
	}
}
