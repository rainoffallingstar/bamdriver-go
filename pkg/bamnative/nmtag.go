package bamnative

import (
	"bytes"
	"fmt"
)

// CalculateNM calculates NM while preserving the historical no-error API.
func CalculateNM(record *Record, reference []byte, isBisulfite bool) int {
	nm, err := CalculateNMChecked(record, reference, isBisulfite)
	if err != nil {
		return 0
	}
	return nm
}

// CalculateNMChecked calculates NM and rejects malformed CIGAR or truncated input.
func CalculateNMChecked(record *Record, reference []byte, isBisulfite bool) (int, error) {
	if record == nil {
		return 0, fmt.Errorf("record is nil")
	}
	if record.RefID < 0 || record.Flags&FlagUnmapped != 0 {
		return 0, nil
	}
	if record.Pos < 0 {
		return 0, fmt.Errorf("record has negative reference position %d", record.Pos)
	}

	readSequence := []byte(record.Seq)
	readIndex := 0
	referenceIndex := int(record.Pos)
	nm := 0

	ensureReadAvailable := func(length int, operation byte) error {
		if length < 0 || readIndex > len(readSequence)-length {
			return fmt.Errorf("CIGAR %d%c exceeds read length at offset %d", length, operation, readIndex)
		}
		return nil
	}
	ensureReferenceAvailable := func(length int, operation byte) error {
		if length < 0 || referenceIndex > len(reference)-length {
			return fmt.Errorf("CIGAR %d%c exceeds reference length at offset %d", length, operation, referenceIndex)
		}
		return nil
	}

	for operationIndex, operation := range record.Cigar {
		if operation.Len <= 0 {
			return 0, fmt.Errorf("CIGAR operation %d has invalid length %d", operationIndex, operation.Len)
		}
		switch operation.Op {
		case CigarMatch:
			if err := ensureReadAvailable(operation.Len, operation.Op); err != nil {
				return 0, err
			}
			if err := ensureReferenceAvailable(operation.Len, operation.Op); err != nil {
				return 0, err
			}
			for baseOffset := 0; baseOffset < operation.Len; baseOffset++ {
				referenceBase := toUpper(reference[referenceIndex+baseOffset])
				readBase := toUpper(readSequence[readIndex+baseOffset])
				isBisulfiteConversion := isBisulfite && ((referenceBase == 'C' && readBase == 'T') || (referenceBase == 'G' && readBase == 'A'))
				if referenceBase != readBase && !isBisulfiteConversion {
					nm++
				}
			}
			readIndex += operation.Len
			referenceIndex += operation.Len
		case CigarEqual:
			if err := ensureReadAvailable(operation.Len, operation.Op); err != nil {
				return 0, err
			}
			if err := ensureReferenceAvailable(operation.Len, operation.Op); err != nil {
				return 0, err
			}
			readIndex += operation.Len
			referenceIndex += operation.Len
		case CigarMismatch:
			if err := ensureReadAvailable(operation.Len, operation.Op); err != nil {
				return 0, err
			}
			if err := ensureReferenceAvailable(operation.Len, operation.Op); err != nil {
				return 0, err
			}
			nm += operation.Len
			readIndex += operation.Len
			referenceIndex += operation.Len
		case CigarInsertion:
			if err := ensureReadAvailable(operation.Len, operation.Op); err != nil {
				return 0, err
			}
			nm += operation.Len
			readIndex += operation.Len
		case CigarDeletion:
			if err := ensureReferenceAvailable(operation.Len, operation.Op); err != nil {
				return 0, err
			}
			nm += operation.Len
			referenceIndex += operation.Len
		case CigarSkip:
			if err := ensureReferenceAvailable(operation.Len, operation.Op); err != nil {
				return 0, err
			}
			referenceIndex += operation.Len
		case CigarSoftClip:
			if err := ensureReadAvailable(operation.Len, operation.Op); err != nil {
				return 0, err
			}
			readIndex += operation.Len
		case CigarHardClip, CigarPadding:
		default:
			return 0, fmt.Errorf("unknown CIGAR operation %q", operation.Op)
		}
	}

	if readIndex != len(readSequence) {
		return 0, fmt.Errorf("CIGAR consumes %d read bases, sequence has %d", readIndex, len(readSequence))
	}
	return nm, nil
}

// toUpper converts byte to uppercase
func toUpper(b byte) byte {
	if b >= 'a' && b <= 'z' {
		return b - 32
	}
	return b
}

// HasNM checks if a record has the NM tag
func HasNM(record *Record, tagName string) bool {
	aux := record.GetAuxField(tagName)
	return aux != nil
}

// FastqToSeq converts a FASTQ-style sequence string to byte slice (for testing)
func FastqToSeq(seq string) []byte {
	return bytes.ToUpper([]byte(seq))
}
