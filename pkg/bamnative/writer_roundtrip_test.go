package bamnative

import (
	"io"
	"os"
	"path/filepath"
	"testing"
)

func TestWriterRoundTripPreservesQualAndAuxArray(t *testing.T) {
	tmp := t.TempDir()
	bamPath := filepath.Join(tmp, "roundtrip.bam")

	h := &Header{
		SortOrder: "coordinate",
		References: []*Reference{
			{ID: 0, Name: "chr1", Len: 1000},
		},
		OtherLines: map[string]string{"@CO": "@CO\troundtrip test"},
	}

	w, err := NewWriter(bamPath, h)
	if err != nil {
		t.Fatalf("NewWriter: %v", err)
	}

	rec := &Record{
		Name:      "read1",
		Flags:     FlagPaired,
		RefID:     0,
		Pos:       10,
		MapQ:      60,
		Cigar:     []CigarOp{{Op: CigarMatch, Len: 5}},
		MateRefID: 0,
		MatePos:   20,
		TLen:      15,
		Seq:       "ACGTN",
		Qual:      nil,
		Aux: []*AuxField{
			{Tag: "XA", Type: AuxTypeArray, ArrayType: AuxTypeUInt8, Value: []byte{1, 2, 3}},
		},
	}

	if err := w.Write(rec); err != nil {
		t.Fatalf("Write: %v", err)
	}
	if err := w.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	f, err := os.Open(bamPath)
	if err != nil {
		t.Fatalf("open BAM: %v", err)
	}
	defer f.Close()

	r, err := NewReader(f)
	if err != nil {
		t.Fatalf("NewReader: %v", err)
	}

	if _, ok := r.Header().OtherLines["@CO"]; !ok {
		t.Fatalf("expected @CO line in header")
	}

	got, err := r.Read()
	if err != nil {
		t.Fatalf("Read: %v", err)
	}

	if got.Seq != rec.Seq {
		t.Fatalf("seq mismatch: got %q want %q", got.Seq, rec.Seq)
	}
	if len(got.Qual) != len(rec.Seq) {
		t.Fatalf("qual length mismatch: got %d want %d", len(got.Qual), len(rec.Seq))
	}
	for i, q := range got.Qual {
		if q != 0xFF {
			t.Fatalf("qual[%d]=%d want 255", i, q)
		}
	}

	aux := got.GetAuxField("XA")
	if aux == nil {
		t.Fatal("missing XA aux field")
	}
	if aux.Type != AuxTypeArray || aux.ArrayType != AuxTypeUInt8 {
		t.Fatalf("unexpected aux type: type=%q arrayType=%q", aux.Type, aux.ArrayType)
	}
	b, ok := aux.Value.([]byte)
	if !ok {
		t.Fatalf("aux value is not []byte: %T", aux.Value)
	}
	if len(b) != 3 || b[0] != 1 || b[1] != 2 || b[2] != 3 {
		t.Fatalf("aux array mismatch: %v", b)
	}

	_, err = r.Read()
	if err != io.EOF {
		t.Fatalf("expected EOF on second read, got %v", err)
	}
}

func TestWriterRejectsQualLengthMismatch(t *testing.T) {
	tmp := t.TempDir()
	bamPath := filepath.Join(tmp, "bad_qual.bam")

	h := &Header{References: []*Reference{{ID: 0, Name: "chr1", Len: 1000}}}
	w, err := NewWriter(bamPath, h)
	if err != nil {
		t.Fatalf("NewWriter: %v", err)
	}
	defer w.Close()

	err = w.Write(&Record{
		Name:  "bad",
		RefID: 0,
		Pos:   1,
		Cigar: []CigarOp{{Op: CigarMatch, Len: 5}},
		Seq:   "ACGTN",
		Qual:  []byte{30, 30},
	})
	if err == nil {
		t.Fatal("expected qual length mismatch error")
	}
}
