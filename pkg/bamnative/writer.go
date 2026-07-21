// Package bamnative provides a pure Go implementation of BAM file parsing
// without relying on external C libraries or biogo/hts.
package bamnative

import (
	"bytes"
	"encoding/binary"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strings"

	"github.com/rainoffallingstar/bamdriver-go/pkg/bgzip"
)

// Writer writes BAM files
type Writer struct {
	header *Header
	bgz    *bgzip.Writer
}

// NewWriter creates a new BAM writer
func NewWriter(path string, header *Header) (*Writer, error) {
	if err := validateHeader(header); err != nil {
		return nil, err
	}

	// Ensure directory exists
	dir := filepath.Dir(path)
	if dir != "." && dir != "" {
		if err := os.MkdirAll(dir, 0755); err != nil {
			return nil, fmt.Errorf("failed to create output directory: %w", err)
		}
	}

	// Create BGZF writer
	bgz, err := bgzip.NewWriter(path)
	if err != nil {
		return nil, fmt.Errorf("failed to create BGZF writer: %w", err)
	}

	w := &Writer{
		header: header,
		bgz:    bgz,
	}

	// Write BAM header
	if err := w.writeHeader(); err != nil {
		bgz.Close()
		return nil, fmt.Errorf("failed to write BAM header: %w", err)
	}

	return w, nil
}

// writeHeader writes the BAM header
func (w *Writer) writeHeader() error {
	// Write BAM magic number first
	magic := [4]byte{'B', 'A', 'M', 1}
	if _, err := w.bgz.Write(magic[:]); err != nil {
		return fmt.Errorf("failed to write BAM magic: %w", err)
	}

	version := strings.TrimSpace(w.header.Version)
	if version == "" {
		version = "1.6"
	}
	sortOrder := strings.TrimSpace(w.header.SortOrder)
	if sortOrder == "" {
		sortOrder = "unknown"
	}

	var headerBuilder strings.Builder
	fmt.Fprintf(&headerBuilder, "@HD\tVN:%s\tSO:%s\n", version, sortOrder)
	for _, reference := range w.header.References {
		fmt.Fprintf(&headerBuilder, "@SQ\tSN:%s\tLN:%d\n", reference.Name, reference.Len)
	}
	if len(w.header.OtherHeaderLines) > 0 {
		for _, headerLine := range w.header.OtherHeaderLines {
			headerLine = strings.TrimSpace(headerLine)
			if headerLine != "" {
				headerBuilder.WriteString(headerLine)
				headerBuilder.WriteByte('\n')
			}
		}
	} else {
		for _, readGroupLine := range w.header.RGLines {
			headerBuilder.WriteString(strings.TrimSpace(readGroupLine))
			headerBuilder.WriteByte('\n')
		}
		for _, programLine := range w.header.PGLines {
			headerBuilder.WriteString(strings.TrimSpace(programLine))
			headerBuilder.WriteByte('\n')
		}
		keys := make([]string, 0, len(w.header.OtherLines))
		for key := range w.header.OtherLines {
			keys = append(keys, key)
		}
		sort.Strings(keys)
		for _, key := range keys {
			headerLine := strings.TrimSpace(w.header.OtherLines[key])
			if headerLine != "" {
				headerBuilder.WriteString(headerLine)
				headerBuilder.WriteByte('\n')
			}
		}
	}
	headerText := headerBuilder.String()

	// Write header length (4 bytes)
	headerLen := int32(len(headerText))
	if err := binary.Write(w.bgz, binary.LittleEndian, headerLen); err != nil {
		return fmt.Errorf("failed to write header length: %w", err)
	}

	// Write header text
	if _, err := w.bgz.Write([]byte(headerText)); err != nil {
		return fmt.Errorf("failed to write header text: %w", err)
	}

	// Write number of reference sequences (4 bytes)
	nRef := int32(len(w.header.References))
	if err := binary.Write(w.bgz, binary.LittleEndian, nRef); err != nil {
		return fmt.Errorf("failed to write reference count: %w", err)
	}

	// Write binary reference data
	for _, ref := range w.header.References {
		// Name length (4 bytes)
		nameLen := int32(len(ref.Name) + 1) // +1 for null terminator
		if err := binary.Write(w.bgz, binary.LittleEndian, nameLen); err != nil {
			return fmt.Errorf("failed to write ref name length: %w", err)
		}

		// Name (null-terminated)
		nameBytes := []byte(ref.Name)
		nameBytes = append(nameBytes, 0)
		if _, err := w.bgz.Write(nameBytes); err != nil {
			return fmt.Errorf("failed to write ref name: %w", err)
		}

		// Reference length (4 bytes)
		if err := binary.Write(w.bgz, binary.LittleEndian, ref.Len); err != nil {
			return fmt.Errorf("failed to write ref length: %w", err)
		}
	}

	return nil
}

// Write validates and encodes a complete BAM record before writing anything to
// the BGZF stream. Validation failures therefore cannot leave a partial record.
func (w *Writer) Write(record *Record) error {
	recordData, err := encodeRecord(record, w.header)
	if err != nil {
		return err
	}
	if len(recordData) > maximumRecordSize {
		return fmt.Errorf("record exceeds maximum size: %d", len(recordData))
	}
	if err := binary.Write(w.bgz, binary.LittleEndian, int32(len(recordData))); err != nil {
		return fmt.Errorf("failed to write block size: %w", err)
	}
	if _, err := w.bgz.Write(recordData); err != nil {
		return fmt.Errorf("failed to write record data: %w", err)
	}
	return nil
}

func encodeRecord(record *Record, header *Header) ([]byte, error) {
	if err := validateRecord(record, header); err != nil {
		return nil, err
	}

	var encoded bytes.Buffer
	queryNameLength := len(record.Name) + 1
	binMQNL := uint32(recordBin(record))<<16 | uint32(record.MapQ)<<8 | uint32(queryNameLength)
	flagNC := uint32(record.Flags)<<16 | uint32(len(record.Cigar))

	coreValues := []interface{}{
		record.RefID,
		record.Pos,
		binMQNL,
		flagNC,
		int32(len(record.Seq)),
		record.MateRefID,
		record.MatePos,
		record.TLen,
	}
	for _, value := range coreValues {
		if err := binary.Write(&encoded, binary.LittleEndian, value); err != nil {
			return nil, fmt.Errorf("encode BAM core: %w", err)
		}
	}

	encoded.WriteString(record.Name)
	encoded.WriteByte(0)
	for _, cigarOperation := range record.Cigar {
		operationCode, err := cigarCharToNum(cigarOperation.Op)
		if err != nil {
			return nil, err
		}
		encodedCigar := uint32(cigarOperation.Len)<<4 | uint32(operationCode)
		if err := binary.Write(&encoded, binary.LittleEndian, encodedCigar); err != nil {
			return nil, fmt.Errorf("encode CIGAR operation: %w", err)
		}
	}
	encoded.Write(encodeSeq(record.Seq))

	if len(record.Qual) == 0 {
		for range record.Seq {
			encoded.WriteByte(0xff)
		}
	} else {
		encoded.Write(record.Qual)
	}

	for _, auxiliaryField := range record.Aux {
		if err := encodeAuxField(&encoded, auxiliaryField); err != nil {
			return nil, fmt.Errorf("encode auxiliary field: %w", err)
		}
	}
	return encoded.Bytes(), nil
}

func validateHeader(header *Header) error {
	if header == nil {
		return fmt.Errorf("header is nil")
	}
	for referenceIndex, reference := range header.References {
		if reference == nil || reference.Name == "" || strings.IndexByte(reference.Name, 0) >= 0 {
			return fmt.Errorf("invalid reference at index %d", referenceIndex)
		}
		if len(reference.Name)+1 > maximumReferenceNameSize || reference.Len < 0 {
			return fmt.Errorf("invalid reference %q", reference.Name)
		}
	}
	return nil
}

func validateRecord(record *Record, header *Header) error {
	if record == nil {
		return fmt.Errorf("record is nil")
	}
	if len(record.Name)+1 > 255 || strings.IndexByte(record.Name, 0) >= 0 {
		return fmt.Errorf("invalid read name length or content")
	}
	if record.RefID < -1 || int(record.RefID) >= len(header.References) {
		return fmt.Errorf("invalid reference ID: %d", record.RefID)
	}
	if record.MateRefID < -1 || int(record.MateRefID) >= len(header.References) {
		return fmt.Errorf("invalid mate reference ID: %d", record.MateRefID)
	}
	if record.RefID >= 0 && record.Flags&FlagUnmapped == 0 && record.Pos < 0 {
		return fmt.Errorf("mapped record has negative position")
	}
	if len(record.Cigar) > 0xffff {
		return fmt.Errorf("too many CIGAR operations: %d", len(record.Cigar))
	}
	readConsumed := 0
	for _, cigarOperation := range record.Cigar {
		if cigarOperation.Len <= 0 || cigarOperation.Len > 0x0fffffff {
			return fmt.Errorf("%w: invalid length %d", ErrInvalidCigar, cigarOperation.Len)
		}
		if _, err := cigarCharToNum(cigarOperation.Op); err != nil {
			return err
		}
		switch cigarOperation.Op {
		case CigarMatch, CigarInsertion, CigarSoftClip, CigarEqual, CigarMismatch:
			readConsumed += cigarOperation.Len
		}
	}
	if len(record.Cigar) > 0 && readConsumed != len(record.Seq) {
		return fmt.Errorf("CIGAR consumes %d read bases, sequence has %d", readConsumed, len(record.Seq))
	}
	if len(record.Qual) != 0 && len(record.Qual) != len(record.Seq) {
		return fmt.Errorf("quality length mismatch: seq=%d qual=%d", len(record.Seq), len(record.Qual))
	}
	for _, base := range []byte(record.Seq) {
		if !isValidSequenceBase(base) {
			return fmt.Errorf("invalid sequence base %q", base)
		}
	}
	for _, auxiliaryField := range record.Aux {
		if err := validateAuxField(auxiliaryField); err != nil {
			return err
		}
	}
	return nil
}

func isValidSequenceBase(base byte) bool {
	const validBases = "=ACMGRSVTWYHKDBNacmgrsvtwyhkdbn"
	return strings.IndexByte(validBases, base) >= 0
}

func validateAuxField(auxiliaryField *AuxField) error {
	if auxiliaryField == nil || len(auxiliaryField.Tag) != 2 {
		return ErrInvalidAux
	}
	for _, character := range []byte(auxiliaryField.Tag) {
		if !((character >= 'A' && character <= 'Z') || (character >= 'a' && character <= 'z') || (character >= '0' && character <= '9')) {
			return ErrInvalidAux
		}
	}
	switch auxiliaryField.Type {
	case AuxTypeChar:
		value, ok := auxiliaryField.Value.(string)
		if !ok || len(value) != 1 {
			return ErrInvalidAux
		}
	case AuxTypeInt8:
		if _, ok := auxiliaryField.Value.(int8); !ok {
			return ErrInvalidAux
		}
	case AuxTypeUInt8:
		if _, ok := auxiliaryField.Value.(uint8); !ok {
			return ErrInvalidAux
		}
	case AuxTypeInt16:
		if _, ok := auxiliaryField.Value.(int16); !ok {
			return ErrInvalidAux
		}
	case AuxTypeUInt16:
		if _, ok := auxiliaryField.Value.(uint16); !ok {
			return ErrInvalidAux
		}
	case AuxTypeInt32:
		if _, ok := auxiliaryField.Value.(int32); !ok {
			return ErrInvalidAux
		}
	case AuxTypeUInt32:
		if _, ok := auxiliaryField.Value.(uint32); !ok {
			return ErrInvalidAux
		}
	case AuxTypeFloat:
		if _, ok := auxiliaryField.Value.(float32); !ok {
			return ErrInvalidAux
		}
	case AuxTypeDouble:
		if _, ok := auxiliaryField.Value.(float64); !ok {
			return ErrInvalidAux
		}
	case AuxTypeString:
		value, ok := auxiliaryField.Value.(string)
		if !ok || strings.IndexByte(value, 0) >= 0 {
			return ErrInvalidAux
		}
	case AuxTypeHex:
		switch value := auxiliaryField.Value.(type) {
		case string:
			if strings.IndexByte(value, 0) >= 0 {
				return ErrInvalidAux
			}
		case []byte:
			if bytes.IndexByte(value, 0) >= 0 {
				return ErrInvalidAux
			}
		default:
			return ErrInvalidAux
		}
	case AuxTypeArray:
		value, ok := auxiliaryField.Value.([]byte)
		elementSize := auxArrayElemSize(auxiliaryField.ArrayType)
		if !ok || elementSize == 0 || len(value)%elementSize != 0 {
			return ErrInvalidAux
		}
	default:
		return ErrInvalidAux
	}
	return nil
}

func encodeAuxField(buffer *bytes.Buffer, auxiliaryField *AuxField) error {
	buffer.WriteString(auxiliaryField.Tag)
	buffer.WriteByte(auxiliaryField.Type)
	switch auxiliaryField.Type {
	case AuxTypeChar:
		buffer.WriteByte(auxiliaryField.Value.(string)[0])
	case AuxTypeInt8:
		buffer.WriteByte(byte(auxiliaryField.Value.(int8)))
	case AuxTypeUInt8:
		buffer.WriteByte(auxiliaryField.Value.(uint8))
	case AuxTypeInt16, AuxTypeUInt16, AuxTypeInt32, AuxTypeUInt32:
		if err := binary.Write(buffer, binary.LittleEndian, auxiliaryField.Value); err != nil {
			return err
		}
	case AuxTypeFloat:
		return binary.Write(buffer, binary.LittleEndian, math.Float32bits(auxiliaryField.Value.(float32)))
	case AuxTypeDouble:
		return binary.Write(buffer, binary.LittleEndian, math.Float64bits(auxiliaryField.Value.(float64)))
	case AuxTypeString:
		buffer.WriteString(auxiliaryField.Value.(string))
		buffer.WriteByte(0)
	case AuxTypeHex:
		switch value := auxiliaryField.Value.(type) {
		case string:
			buffer.WriteString(value)
		case []byte:
			buffer.Write(value)
		}
		buffer.WriteByte(0)
	case AuxTypeArray:
		value := auxiliaryField.Value.([]byte)
		elementSize := auxArrayElemSize(auxiliaryField.ArrayType)
		buffer.WriteByte(auxiliaryField.ArrayType)
		if err := binary.Write(buffer, binary.LittleEndian, int32(len(value)/elementSize)); err != nil {
			return err
		}
		buffer.Write(value)
	}
	return nil
}

// encodeSeq encodes a DNA sequence to 4-bit BAM format
func encodeSeq(seq string) []byte {
	seqLen := len(seq)
	if seqLen == 0 {
		return []byte{}
	}

	seqBytes := (seqLen + 1) / 2
	result := make([]byte, seqBytes)

	for i := 0; i < seqLen; i++ {
		base := seqCharToNum(seq[i])
		if i%2 == 0 {
			result[i/2] = base << 4
		} else {
			result[i/2] |= base
		}
	}

	return result
}

func seqCharToNum(base byte) byte {
	switch base {
	case '=', 0:
		return 0
	case 'A', 'a':
		return 1
	case 'C', 'c':
		return 2
	case 'M', 'm':
		return 3
	case 'G', 'g':
		return 4
	case 'R', 'r':
		return 5
	case 'S', 's':
		return 6
	case 'V', 'v':
		return 7
	case 'T', 't':
		return 8
	case 'W', 'w':
		return 9
	case 'Y', 'y':
		return 10
	case 'H', 'h':
		return 11
	case 'K', 'k':
		return 12
	case 'D', 'd':
		return 13
	case 'B', 'b':
		return 14
	case 'N', 'n':
		return 15
	default:
		return 15
	}
}

func cigarCharToNum(operation byte) (byte, error) {
	switch operation {
	case CigarMatch:
		return 0, nil
	case CigarInsertion:
		return 1, nil
	case CigarDeletion:
		return 2, nil
	case CigarSkip:
		return 3, nil
	case CigarSoftClip:
		return 4, nil
	case CigarHardClip:
		return 5, nil
	case CigarPadding:
		return 6, nil
	case CigarEqual:
		return 7, nil
	case CigarMismatch:
		return 8, nil
	default:
		return 0, fmt.Errorf("%w: %q", ErrInvalidCigar, operation)
	}
}

func recordBin(record *Record) uint16 {
	if record.RefID < 0 || len(record.Cigar) == 0 {
		return 0
	}
	beg := int(record.Pos)
	end := beg
	for _, op := range record.Cigar {
		switch op.Op {
		case CigarMatch, CigarDeletion, CigarSkip, CigarEqual, CigarMismatch:
			end += op.Len
		}
	}
	if end <= beg {
		end = beg + 1
	}
	return uint16(reg2bin(beg, end))
}

// Close closes the writer
func (w *Writer) Close() error {
	return w.bgz.Close()
}

// VirtualOffset returns the BAI virtual offset of the next byte to be written
// to the underlying BGZF stream.
func (w *Writer) VirtualOffset() int64 {
	return w.bgz.VirtualOffset()
}

// WriterAt is a wrapper around Writer that allows writing at a specific position
type WriterAt struct {
	*Writer
}

// WriteRecords writes multiple records
func (w *Writer) WriteRecords(records []*Record) error {
	for _, record := range records {
		if err := w.Write(record); err != nil {
			return err
		}
	}
	return nil
}
