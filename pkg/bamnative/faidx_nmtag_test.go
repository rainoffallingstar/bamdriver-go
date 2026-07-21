package bamnative

import (
	"os"
	"path/filepath"
	"sync"
	"testing"
)

func TestFastaIndexUsesPhysicalLineWidths(t *testing.T) {
	testCases := []struct {
		name         string
		content      string
		wantIndex    string
		wantSequence map[string]string
	}{
		{
			name:      "LF",
			content:   ">chr1 description\nACGT\nTG\n>chr2\nNNNN",
			wantIndex: "chr1\t6\t18\t4\t5\nchr2\t4\t32\t4\t4\n",
			wantSequence: map[string]string{
				"chr1": "ACGTTG",
				"chr2": "NNNN",
			},
		},
		{
			name:      "CRLF",
			content:   ">chr1\r\nACGT\r\nTG\r\n>chr2\r\nNNNN",
			wantIndex: "chr1\t6\t7\t4\t6\nchr2\t4\t24\t4\t4\n",
			wantSequence: map[string]string{
				"chr1": "ACGTTG",
				"chr2": "NNNN",
			},
		},
		{
			name:      "full final line without newline",
			content:   ">chr1\nACGT\nTGCA",
			wantIndex: "chr1\t8\t6\t4\t5\n",
			wantSequence: map[string]string{
				"chr1": "ACGTTGCA",
			},
		},
	}

	for _, testCase := range testCases {
		t.Run(testCase.name, func(t *testing.T) {
			temporaryDirectory := t.TempDir()
			fastaPath := filepath.Join(temporaryDirectory, "reference.fa")
			if err := os.WriteFile(fastaPath, []byte(testCase.content), 0o600); err != nil {
				t.Fatalf("WriteFile FASTA: %v", err)
			}

			reader, err := NewFastaReader(fastaPath)
			if err != nil {
				t.Fatalf("NewFastaReader: %v", err)
			}
			indexData, err := os.ReadFile(fastaPath + ".fai")
			if err != nil {
				t.Fatalf("ReadFile index: %v", err)
			}
			if string(indexData) != testCase.wantIndex {
				t.Fatalf("index = %q, want %q", indexData, testCase.wantIndex)
			}
			for name, wantSequence := range testCase.wantSequence {
				sequence, ok := reader.GetSequence(name)
				if !ok {
					t.Fatalf("GetSequence(%q) returned not found", name)
				}
				if string(sequence) != wantSequence {
					t.Fatalf("GetSequence(%q) = %q, want %q", name, sequence, wantSequence)
				}
			}
		})
	}
}

func TestFastaReaderRebuildsInvalidExistingIndex(t *testing.T) {
	temporaryDirectory := t.TempDir()
	fastaPath := filepath.Join(temporaryDirectory, "reference.fa")
	if err := os.WriteFile(fastaPath, []byte(">chr1\nACGT\nTG\n"), 0o600); err != nil {
		t.Fatalf("WriteFile FASTA: %v", err)
	}
	if err := os.WriteFile(fastaPath+".fai", []byte("chr1\t6\t6\t0\t5\n"), 0o600); err != nil {
		t.Fatalf("WriteFile invalid index: %v", err)
	}

	reader, err := NewFastaReader(fastaPath)
	if err != nil {
		t.Fatalf("NewFastaReader: %v", err)
	}
	sequence, ok := reader.GetSequence("chr1")
	if !ok || string(sequence) != "ACGTTG" {
		t.Fatalf("GetSequence = %q, %v", sequence, ok)
	}
	indexData, err := os.ReadFile(fastaPath + ".fai")
	if err != nil {
		t.Fatalf("ReadFile rebuilt index: %v", err)
	}
	if string(indexData) != "chr1\t6\t6\t4\t5\n" {
		t.Fatalf("rebuilt index = %q", indexData)
	}
}

func TestFastaReaderConcurrentGetSequence(t *testing.T) {
	temporaryDirectory := t.TempDir()
	fastaPath := filepath.Join(temporaryDirectory, "reference.fa")
	if err := os.WriteFile(fastaPath, []byte(">chr1\nACGT\n>chr2\nTGCA\n"), 0o600); err != nil {
		t.Fatalf("WriteFile FASTA: %v", err)
	}
	reader, err := NewFastaReader(fastaPath)
	if err != nil {
		t.Fatalf("NewFastaReader: %v", err)
	}

	var waitGroup sync.WaitGroup
	for workerIndex := 0; workerIndex < 16; workerIndex++ {
		waitGroup.Add(1)
		go func() {
			defer waitGroup.Done()
			for iteration := 0; iteration < 100; iteration++ {
				sequence, ok := reader.GetSequence("chr1")
				if !ok || string(sequence) != "ACGT" {
					t.Errorf("GetSequence = %q, %v", sequence, ok)
					return
				}
			}
		}()
	}
	waitGroup.Wait()
}

func TestCalculateNMChecked(t *testing.T) {
	record := &Record{
		RefID: 0,
		Pos:   0,
		Seq:   "TACGGA",
		Cigar: []CigarOp{
			{Op: CigarSoftClip, Len: 1},
			{Op: CigarMatch, Len: 3},
			{Op: CigarInsertion, Len: 1},
			{Op: CigarMismatch, Len: 1},
			{Op: CigarDeletion, Len: 1},
		},
	}

	nm, err := CalculateNMChecked(record, []byte("ACGTA"), false)
	if err != nil {
		t.Fatalf("CalculateNMChecked: %v", err)
	}
	if nm != 3 {
		t.Fatalf("NM = %d, want 3", nm)
	}
}

func TestCalculateNMCheckedRejectsTruncatedInputs(t *testing.T) {
	testCases := []struct {
		name      string
		record    *Record
		reference []byte
	}{
		{
			name: "read too short",
			record: &Record{
				RefID: 0,
				Pos:   0,
				Seq:   "AC",
				Cigar: []CigarOp{{Op: CigarMatch, Len: 3}},
			},
			reference: []byte("ACG"),
		},
		{
			name: "reference too short",
			record: &Record{
				RefID: 0,
				Pos:   2,
				Seq:   "AC",
				Cigar: []CigarOp{{Op: CigarMatch, Len: 2}},
			},
			reference: []byte("ACG"),
		},
		{
			name: "unknown operation",
			record: &Record{
				RefID: 0,
				Pos:   0,
				Seq:   "A",
				Cigar: []CigarOp{{Op: 'Q', Len: 1}},
			},
			reference: []byte("A"),
		},
		{
			name: "zero length operation",
			record: &Record{
				RefID: 0,
				Pos:   0,
				Seq:   "A",
				Cigar: []CigarOp{{Op: CigarMatch, Len: 0}},
			},
			reference: []byte("A"),
		},
	}

	for _, testCase := range testCases {
		t.Run(testCase.name, func(t *testing.T) {
			if _, err := CalculateNMChecked(testCase.record, testCase.reference, false); err == nil {
				t.Fatal("CalculateNMChecked succeeded for malformed input")
			}
		})
	}
}

func TestCalculateNMCheckedBisulfiteConversions(t *testing.T) {
	record := &Record{
		RefID: 0,
		Pos:   0,
		Seq:   "TA",
		Cigar: []CigarOp{{Op: CigarMatch, Len: 2}},
	}
	withoutBisulfite, err := CalculateNMChecked(record, []byte("CG"), false)
	if err != nil {
		t.Fatalf("CalculateNMChecked standard: %v", err)
	}
	withBisulfite, err := CalculateNMChecked(record, []byte("CG"), true)
	if err != nil {
		t.Fatalf("CalculateNMChecked bisulfite: %v", err)
	}
	if withoutBisulfite != 2 || withBisulfite != 0 {
		t.Fatalf("NM standard/bisulfite = %d/%d, want 2/0", withoutBisulfite, withBisulfite)
	}
}
