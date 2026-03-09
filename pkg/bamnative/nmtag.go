package bamnative

import (
	"bytes"
)

// CalculateNM calculates the NM (edit distance) tag by comparing read to reference
// NM = mismatches + insertions + deletions
func CalculateNM(record *Record, ref []byte, isBisulfite bool) int {
	if len(record.Seq) == 0 || record.RefID < 0 {
		return 0
	}

	// record.Seq is already decoded ASCII by the BAM reader; use it directly
	readSeq := []byte(record.Seq)

	// Get reference position (0-based in BAM)
	refPos := int(record.Pos)

	// Parse CIGAR and calculate NM
	nm := 0

	// Position in read sequence
	readIdx := 0
	// Position in reference sequence
	refIdx := refPos

	for _, cigar := range record.Cigar {
		switch cigar.Op {
		case CigarMatch:
			// M: match or mismatch - count as full length
			// Need to compare each position to find actual mismatches
			for i := 0; i < cigar.Len; i++ {
				if refIdx+i >= len(ref) || readIdx+i >= len(readSeq) {
					break
				}
				refBase := toUpper(ref[refIdx+i])
				readBase := toUpper(readSeq[readIdx+i])

				if isBisulfite {
					// In bisulfite mode, C->T on read is not counted as mismatch
					// if reference has C
					if refBase == 'C' && readBase == 'T' {
						// forward strand: C->T conversion, not a mismatch
					} else if refBase == 'G' && readBase == 'A' {
						// reverse strand: G->A conversion (complement of C->T), not a mismatch
					} else if refBase != readBase {
						nm++
					}
				} else {
					if refBase != readBase {
						nm++
					}
				}
			}
			readIdx += cigar.Len
			refIdx += cigar.Len

		case CigarEqual:
			// =: exact match - no mismatches
			readIdx += cigar.Len
			refIdx += cigar.Len

		case CigarMismatch:
			// X: all are mismatches
			nm += cigar.Len
			readIdx += cigar.Len
			refIdx += cigar.Len

		case CigarInsertion:
			// I: insertion relative to reference
			nm += cigar.Len
			readIdx += cigar.Len

		case CigarDeletion:
			// D: deletion relative to reference
			nm += cigar.Len
			refIdx += cigar.Len

		case CigarSkip:
			// N: intron skip - not counted in NM per SAM spec
			refIdx += cigar.Len

		case CigarSoftClip:
			// S: soft clip - not counted in NM; tracked separately as clips
			readIdx += cigar.Len

		case CigarHardClip:
			// H: hard clip - sequence not present in read
			// Does not affect NM - no position changes

		case CigarPadding:
			// P: padding - no effect on NM
		}
	}

	return nm
}

// decodeSeq converts BAM 4-bit encoded sequence to byte slice
// BAM uses: 0=?, 1=A, 2=C, 3=G, 4=T, 5=N, 6=<, 7=>
func decodeSeq(bamSeq string) []byte {
	// BAM sequence is stored as a string where each character encodes 2 bases
	// using 4-bit encoding: 0-15 per nibble
	seqMap := map[byte]byte{
		'0': '=', // or '?'
		'1': 'A',
		'2': 'C',
		'3': 'G',
		'4': 'T',
		'5': 'N',
		'6': '<',
		'7': '>',
		'8': '?',
		'9': '?',
		'a': 'A', // Could be lower case?
		'b': 'C',
		'c': 'G',
		'd': 'T',
		'e': 'N',
		'f': '?',
	}

	result := make([]byte, 0, len(bamSeq)*2)

	for i := 0; i < len(bamSeq); i++ {
		// Each byte encodes two bases
		hi := (bamSeq[i] >> 4) & 0x0F
		lo := bamSeq[i] & 0x0F

		if hi <= 9 {
			if hi <= 5 {
				result = append(result, seqMap[byte('0'+hi)])
			} else {
				result = append(result, seqMap[byte('a'+hi-10)])
			}
		}

		if lo <= 9 {
			if lo <= 5 {
				result = append(result, seqMap[byte('0'+lo)])
			} else {
				result = append(result, seqMap[byte('a'+lo-10)])
			}
		}
	}

	return result
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
