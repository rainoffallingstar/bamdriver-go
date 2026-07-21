package bamnative

import (
	"bytes"
	"io"
	"testing"
)

func FuzzBAMReaderMalformedInput(fuzz *testing.F) {
	fuzz.Add([]byte{})
	fuzz.Add([]byte("BAM\x01"))

	fuzz.Fuzz(func(t *testing.T, bamData []byte) {
		reader, err := NewReader(bytes.NewReader(bamData))
		if err != nil {
			return
		}
		for recordCount := 0; recordCount < 1024; recordCount++ {
			_, readErr := reader.Read()
			if readErr != nil {
				if readErr == io.EOF {
					return
				}
				return
			}
		}
		t.Fatal("reader accepted more than 1024 records from bounded fuzz input")
	})
}

func FuzzCalculateNMChecked(fuzz *testing.F) {
	fuzz.Add("ACGT", []byte("ACGT"), uint8(CigarMatch), uint16(4), int32(0))
	fuzz.Fuzz(func(t *testing.T, sequence string, reference []byte, operation byte, operationLength uint16, position int32) {
		if len(sequence) > 4096 || len(reference) > 4096 {
			t.Skip()
		}
		record := &Record{
			RefID: 0,
			Pos:   position,
			Seq:   sequence,
			Cigar: []CigarOp{{Op: operation, Len: int(operationLength)}},
		}
		_, _ = CalculateNMChecked(record, reference, false)
	})
}
