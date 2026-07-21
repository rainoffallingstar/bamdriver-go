package bamnative

import (
	"container/heap"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"sort"
)

// SortOptions contains options for BAM sorting
type SortOptions struct {
	OutputPath         string
	ByName             bool
	MemoryLimitBytes   int64
	TemporaryDirectory string
}

const (
	defaultSortMemoryLimit = int64(64 << 20)
	maximumMergeFanIn      = 64
)

// Sort performs an external stable sort bounded by the configured memory limit.
func Sort(inputPath string, options *SortOptions) error {
	if options == nil || options.OutputPath == "" {
		return fmt.Errorf("sort output path is required")
	}
	memoryLimit := options.MemoryLimitBytes
	if memoryLimit <= 0 {
		memoryLimit = defaultSortMemoryLimit
	}
	temporaryRoot := options.TemporaryDirectory
	if temporaryRoot == "" {
		temporaryRoot = filepath.Dir(options.OutputPath)
	}
	temporaryDirectory, err := os.MkdirTemp(temporaryRoot, ".bam-sort-*")
	if err != nil {
		return fmt.Errorf("failed to create sort temporary directory: %w", err)
	}
	defer os.RemoveAll(temporaryDirectory)

	inputFile, err := os.Open(inputPath)
	if err != nil {
		return fmt.Errorf("failed to open input file: %w", err)
	}
	defer inputFile.Close()
	reader, err := NewReader(inputFile)
	if err != nil {
		return fmt.Errorf("failed to create reader: %w", err)
	}
	header := reader.Header()
	if options.ByName {
		header.SortOrder = "queryname"
	} else {
		header.SortOrder = "coordinate"
	}

	runPaths, err := createSortedRuns(reader, header, options.ByName, memoryLimit, temporaryDirectory)
	if err != nil {
		return err
	}
	if len(runPaths) == 0 {
		writer, err := NewWriter(options.OutputPath, header)
		if err != nil {
			return err
		}
		return writer.Close()
	}

	mergePass := 0
	for len(runPaths) > maximumMergeFanIn {
		var mergedRunPaths []string
		for groupStart := 0; groupStart < len(runPaths); groupStart += maximumMergeFanIn {
			groupEnd := groupStart + maximumMergeFanIn
			if groupEnd > len(runPaths) {
				groupEnd = len(runPaths)
			}
			mergedPath := filepath.Join(
				temporaryDirectory,
				fmt.Sprintf("merge-%03d-%06d.bam", mergePass, len(mergedRunPaths)),
			)
			if err := mergeSortedRuns(runPaths[groupStart:groupEnd], mergedPath, header, options.ByName); err != nil {
				return err
			}
			mergedRunPaths = append(mergedRunPaths, mergedPath)
		}
		for _, runPath := range runPaths {
			_ = os.Remove(runPath)
		}
		runPaths = mergedRunPaths
		mergePass++
	}

	return mergeSortedRuns(runPaths, options.OutputPath, header, options.ByName)
}

func createSortedRuns(
	reader *Reader,
	header *Header,
	byName bool,
	memoryLimit int64,
	temporaryDirectory string,
) ([]string, error) {
	var runPaths []string
	var records []*Record
	var estimatedBytes int64

	flushRun := func() error {
		if len(records) == 0 {
			return nil
		}
		sort.SliceStable(records, func(leftIndex, rightIndex int) bool {
			return compareSortRecords(records[leftIndex], records[rightIndex], byName) < 0
		})
		runPath := filepath.Join(temporaryDirectory, fmt.Sprintf("run-%06d.bam", len(runPaths)))
		writer, err := NewWriter(runPath, header)
		if err != nil {
			return fmt.Errorf("failed to create sorted run: %w", err)
		}
		for _, record := range records {
			if err := writer.Write(record); err != nil {
				_ = writer.Close()
				return fmt.Errorf("failed to write sorted run: %w", err)
			}
		}
		if err := writer.Close(); err != nil {
			return fmt.Errorf("failed to close sorted run: %w", err)
		}
		runPaths = append(runPaths, runPath)
		records = nil
		estimatedBytes = 0
		return nil
	}

	for {
		record, readErr := reader.Read()
		if readErr != nil {
			if readErr == io.EOF {
				break
			}
			return nil, fmt.Errorf("failed to read record while sorting: %w", readErr)
		}
		if err := validateRecord(record, header); err != nil {
			return nil, fmt.Errorf("invalid record while sorting: %w", err)
		}
		recordBytes := estimateRecordMemory(record)
		if len(records) > 0 && estimatedBytes+recordBytes > memoryLimit {
			if err := flushRun(); err != nil {
				return nil, err
			}
		}
		records = append(records, record)
		estimatedBytes += recordBytes
	}
	if err := flushRun(); err != nil {
		return nil, err
	}
	return runPaths, nil
}

func mergeSortedRuns(runPaths []string, outputPath string, header *Header, byName bool) error {
	openRuns := make([]*mergeRun, 0, len(runPaths))
	defer func() {
		for _, run := range openRuns {
			_ = run.file.Close()
		}
	}()

	queue := &mergeRunHeap{byName: byName}
	for sourceIndex, runPath := range runPaths {
		runFile, err := os.Open(runPath)
		if err != nil {
			return fmt.Errorf("failed to open sorted run: %w", err)
		}
		runReader, err := NewReader(runFile)
		if err != nil {
			_ = runFile.Close()
			return fmt.Errorf("failed to read sorted run header: %w", err)
		}
		run := &mergeRun{file: runFile, reader: runReader, sourceIndex: sourceIndex}
		openRuns = append(openRuns, run)
		record, readErr := runReader.Read()
		if readErr == io.EOF {
			continue
		}
		if readErr != nil {
			return fmt.Errorf("failed to read sorted run: %w", readErr)
		}
		run.record = record
		heap.Push(queue, run)
	}

	writer, err := NewWriter(outputPath, header)
	if err != nil {
		return fmt.Errorf("failed to create merged BAM: %w", err)
	}
	for queue.Len() > 0 {
		run := heap.Pop(queue).(*mergeRun)
		if err := writer.Write(run.record); err != nil {
			_ = writer.Close()
			return fmt.Errorf("failed to write merged BAM: %w", err)
		}
		nextRecord, readErr := run.reader.Read()
		if readErr == nil {
			run.record = nextRecord
			heap.Push(queue, run)
			continue
		}
		if readErr != io.EOF {
			_ = writer.Close()
			return fmt.Errorf("failed to advance sorted run: %w", readErr)
		}
	}
	if err := writer.Close(); err != nil {
		return fmt.Errorf("failed to close merged BAM: %w", err)
	}
	return nil
}

type mergeRun struct {
	file        *os.File
	reader      *Reader
	record      *Record
	sourceIndex int
}

type mergeRunHeap struct {
	items  []*mergeRun
	byName bool
}

func (queue mergeRunHeap) Len() int { return len(queue.items) }

func (queue mergeRunHeap) Less(leftIndex, rightIndex int) bool {
	comparison := compareSortRecords(
		queue.items[leftIndex].record,
		queue.items[rightIndex].record,
		queue.byName,
	)
	if comparison == 0 {
		return queue.items[leftIndex].sourceIndex < queue.items[rightIndex].sourceIndex
	}
	return comparison < 0
}

func (queue mergeRunHeap) Swap(leftIndex, rightIndex int) {
	queue.items[leftIndex], queue.items[rightIndex] = queue.items[rightIndex], queue.items[leftIndex]
}

func (queue *mergeRunHeap) Push(value interface{}) {
	queue.items = append(queue.items, value.(*mergeRun))
}

func (queue *mergeRunHeap) Pop() interface{} {
	lastIndex := len(queue.items) - 1
	value := queue.items[lastIndex]
	queue.items[lastIndex] = nil
	queue.items = queue.items[:lastIndex]
	return value
}

func compareSortRecords(leftRecord, rightRecord *Record, byName bool) int {
	if byName {
		if leftRecord.Name < rightRecord.Name {
			return -1
		}
		if leftRecord.Name > rightRecord.Name {
			return 1
		}
		return 0
	}
	return compareCoordinateRecords(leftRecord, rightRecord)
}

func estimateRecordMemory(record *Record) int64 {
	estimatedBytes := int64(192 + len(record.Name) + len(record.Seq) + len(record.Qual))
	estimatedBytes += int64(len(record.Cigar) * 16)
	estimatedBytes += int64(len(record.Aux) * 48)
	for _, auxiliaryField := range record.Aux {
		if auxiliaryField == nil {
			continue
		}
		estimatedBytes += int64(len(auxiliaryField.Tag))
		switch value := auxiliaryField.Value.(type) {
		case string:
			estimatedBytes += int64(len(value))
		case []byte:
			estimatedBytes += int64(len(value))
		}
	}
	return estimatedBytes
}

// IsSorted checks the records rather than trusting the header sort-order field.
func IsSorted(path string) (bool, error) {
	inputFile, err := os.Open(path)
	if err != nil {
		return false, fmt.Errorf("failed to open file: %w", err)
	}
	defer inputFile.Close()

	reader, err := NewReader(inputFile)
	if err != nil {
		return false, fmt.Errorf("failed to create reader: %w", err)
	}

	referenceCount := len(reader.Header().References)
	var previousRecord *Record
	for {
		record, readErr := reader.Read()
		if readErr != nil {
			if readErr == io.EOF {
				return true, nil
			}
			return false, fmt.Errorf("failed to read record: %w", readErr)
		}
		if err := validateSortCoordinates(record, referenceCount); err != nil {
			return false, err
		}
		if previousRecord != nil && compareCoordinateRecords(previousRecord, record) > 0 {
			return false, nil
		}
		previousRecord = record
	}
}

// SortAndIndexIfNeeded produces a coordinate-sorted BAM and its BAI index.
func SortAndIndexIfNeeded(inputPath, outputPath string) error {
	isSorted, err := IsSorted(inputPath)
	if err != nil {
		return err
	}

	if isSorted {
		sameFile, err := pathsReferToSameFile(inputPath, outputPath)
		if err != nil {
			return err
		}
		if !sameFile {
			if err := copyBamFile(inputPath, outputPath); err != nil {
				return fmt.Errorf("failed to copy sorted BAM: %w", err)
			}
		}
	} else {
		if err := Sort(inputPath, &SortOptions{OutputPath: outputPath}); err != nil {
			return err
		}
	}

	if err := BuildIndex(outputPath); err != nil {
		return fmt.Errorf("failed to build BAM index: %w", err)
	}
	return nil
}

func compareCoordinateRecords(leftRecord, rightRecord *Record) int {
	leftHasCoordinate := recordHasCoordinate(leftRecord)
	rightHasCoordinate := recordHasCoordinate(rightRecord)
	if leftHasCoordinate != rightHasCoordinate {
		if leftHasCoordinate {
			return -1
		}
		return 1
	}
	if !leftHasCoordinate {
		return 0
	}
	if leftRecord.RefID < rightRecord.RefID {
		return -1
	}
	if leftRecord.RefID > rightRecord.RefID {
		return 1
	}
	if leftRecord.Pos < rightRecord.Pos {
		return -1
	}
	if leftRecord.Pos > rightRecord.Pos {
		return 1
	}
	return 0
}

func recordHasCoordinate(record *Record) bool {
	return record.RefID >= 0 && record.Pos >= 0
}

func validateSortCoordinates(record *Record, referenceCount int) error {
	if record.RefID < -1 || int(record.RefID) >= referenceCount {
		return fmt.Errorf("record %q has invalid reference ID %d", record.Name, record.RefID)
	}
	if record.Pos < -1 {
		return fmt.Errorf("record %q has invalid position %d", record.Name, record.Pos)
	}
	if record.Flags&FlagUnmapped == 0 && !recordHasCoordinate(record) {
		return fmt.Errorf("mapped record %q has no coordinate", record.Name)
	}
	return nil
}

func pathsReferToSameFile(firstPath, secondPath string) (bool, error) {
	firstInfo, err := os.Stat(firstPath)
	if err != nil {
		return false, fmt.Errorf("failed to inspect input BAM: %w", err)
	}
	secondInfo, secondErr := os.Stat(secondPath)
	if secondErr == nil {
		return os.SameFile(firstInfo, secondInfo), nil
	}
	if !os.IsNotExist(secondErr) {
		return false, fmt.Errorf("failed to inspect output BAM: %w", secondErr)
	}
	firstAbsolute, err := filepath.Abs(firstPath)
	if err != nil {
		return false, fmt.Errorf("failed to resolve input BAM path: %w", err)
	}
	secondAbsolute, err := filepath.Abs(secondPath)
	if err != nil {
		return false, fmt.Errorf("failed to resolve output BAM path: %w", err)
	}
	return filepath.Clean(firstAbsolute) == filepath.Clean(secondAbsolute), nil
}

func copyBamFile(sourcePath, destinationPath string) error {
	sourceFile, err := os.Open(sourcePath)
	if err != nil {
		return err
	}
	defer sourceFile.Close()

	destinationFile, err := os.Create(destinationPath)
	if err != nil {
		return err
	}
	if _, err := io.Copy(destinationFile, sourceFile); err != nil {
		_ = destinationFile.Close()
		return err
	}
	return destinationFile.Close()
}
