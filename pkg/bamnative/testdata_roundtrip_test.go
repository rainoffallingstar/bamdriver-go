package bamnative

import (
	"bytes"
	"fmt"
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"reflect"
	"testing"
)

func TestRoundTripTestdataBAMPreservesRecords(t *testing.T) {
	t.Parallel()

	for _, name := range []string{
		"Test_hg19_NRAS.bam",
		"Test_mm10_NRAS.bam",
	} {
		t.Run(name, func(t *testing.T) {
			t.Parallel()

			inPath := filepath.Join("..", "..", "testdata", name)
			inF, err := os.Open(inPath)
			if err != nil {
				t.Fatalf("open input BAM: %v", err)
			}
			defer inF.Close()

			inR, err := NewReader(inF)
			if err != nil {
				t.Fatalf("NewReader(input): %v", err)
			}

			hdr := inR.Header()
			tmp := t.TempDir()
			outPath := filepath.Join(tmp, "roundtrip.bam")

			w, err := NewWriter(outPath, hdr)
			if err != nil {
				t.Fatalf("NewWriter: %v", err)
			}

			var want []*Record
			for {
				rec, err := inR.Read()
				if err == io.EOF {
					break
				}
				if err != nil {
					_ = w.Close()
					t.Fatalf("Read(input): %v", err)
				}
				want = append(want, rec)
				if err := w.Write(rec); err != nil {
					_ = w.Close()
					t.Fatalf("Write(output): %v", err)
				}
			}
			if err := w.Close(); err != nil {
				t.Fatalf("Close(output): %v", err)
			}

			outF, err := os.Open(outPath)
			if err != nil {
				t.Fatalf("open output BAM: %v", err)
			}
			defer outF.Close()

			outR, err := NewReader(outF)
			if err != nil {
				t.Fatalf("NewReader(output): %v", err)
			}

			var got []*Record
			for {
				rec, err := outR.Read()
				if err == io.EOF {
					break
				}
				if err != nil {
					t.Fatalf("Read(output): %v", err)
				}
				got = append(got, rec)
			}

			if len(got) != len(want) {
				t.Fatalf("record count mismatch: got %d want %d", len(got), len(want))
			}
			for i := range want {
				if err := compareRecords(got[i], want[i]); err != nil {
					t.Fatalf("record[%d] mismatch: %v", i, err)
				}
			}
		})
	}
}

func TestWriterOutputIsSamtoolsCompatible(t *testing.T) {
	t.Parallel()

	if _, err := exec.LookPath("samtools"); err != nil {
		t.Skip("samtools not available")
	}

	tmp := t.TempDir()
	outPath := filepath.Join(tmp, "samtools_compat.bam")

	h := &Header{
		SortOrder: "coordinate",
		References: []*Reference{
			{ID: 0, Name: "chr1", Len: 1000},
		},
		OtherLines: map[string]string{"@CO": "@CO\tsamtools compat test"},
	}

	w, err := NewWriter(outPath, h)
	if err != nil {
		t.Fatalf("NewWriter: %v", err)
	}

	for i := 0; i < 5; i++ {
		rec := &Record{
			Name:      "read",
			Flags:     0,
			RefID:     0,
			Pos:       int32(10 + i),
			MapQ:      60,
			Cigar:     []CigarOp{{Op: CigarMatch, Len: 5}},
			MateRefID: -1,
			MatePos:   -1,
			TLen:      0,
			Seq:       "ACGTN",
			Qual:      []byte{30, 31, 32, 33, 34},
		}
		if err := w.Write(rec); err != nil {
			_ = w.Close()
			t.Fatalf("Write: %v", err)
		}
	}
	if err := w.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	runSamtools(t, "quickcheck", outPath)
	runSamtools(t, "view", "-H", outPath)
	runSamtools(t, "view", "-c", outPath)
}

func runSamtools(t *testing.T, args ...string) {
	t.Helper()
	cmd := exec.Command("samtools", args...)
	var stderr bytes.Buffer
	cmd.Stderr = &stderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("samtools %v failed: %v; stderr: %s", args, err, stderr.String())
	}
}

func compareRecords(got, want *Record) error {
	if got.Name != want.Name {
		return fmt.Errorf("name got %q want %q", got.Name, want.Name)
	}
	if got.Flags != want.Flags {
		return fmt.Errorf("flags got %d want %d", got.Flags, want.Flags)
	}
	if got.RefID != want.RefID || got.Pos != want.Pos {
		return fmt.Errorf("pos got (%d,%d) want (%d,%d)", got.RefID, got.Pos, want.RefID, want.Pos)
	}
	if got.MapQ != want.MapQ {
		return fmt.Errorf("mapq got %d want %d", got.MapQ, want.MapQ)
	}
	if got.MateRefID != want.MateRefID || got.MatePos != want.MatePos {
		return fmt.Errorf("mate got (%d,%d) want (%d,%d)", got.MateRefID, got.MatePos, want.MateRefID, want.MatePos)
	}
	if got.TLen != want.TLen {
		return fmt.Errorf("tlen got %d want %d", got.TLen, want.TLen)
	}
	if got.Seq != want.Seq {
		return fmt.Errorf("seq got %q want %q", got.Seq, want.Seq)
	}
	if !bytes.Equal(got.Qual, want.Qual) {
		return fmt.Errorf("qual mismatch")
	}
	if !reflect.DeepEqual(got.Cigar, want.Cigar) {
		return fmt.Errorf("cigar mismatch got %v want %v", got.Cigar, want.Cigar)
	}
	if !reflect.DeepEqual(got.Aux, want.Aux) {
		return fmt.Errorf("aux mismatch")
	}
	return nil
}
