// Package bamnative provides a pure Go implementation of BAM file parsing
// without relying on external C libraries or biogo/hts.
package bamnative

import (
	"encoding/binary"
	"errors"
	"fmt"
	"io"
	"math"
	"strings"

	"github.com/rainoffallingstar/bamdriver-go/pkg/bgzip"
)

var (
	ErrInvalidHeader     = errors.New("invalid BAM header")
	ErrInvalidRecord     = errors.New("invalid BAM record")
	ErrMissingField      = errors.New("missing required field")
	ErrInvalidCigar      = errors.New("invalid CIGAR operation")
	ErrInvalidAux        = errors.New("invalid auxiliary field")
	ErrNoRecords         = errors.New("no records in BAM")
	ErrReferenceNotFound = errors.New("reference not found")
)

// Header represents the BAM file header
type Header struct {
	Version    string            // Version string (e.g., "1.0" or empty)
	SortOrder  string            // Sort order (e.g., "coordinate", "queryname", "unknown")
	OtherLines map[string]string // Other header lines (non-@RG/@PG lines)
	References []*Reference      // Reference sequences
	RGLines    []string          // @RG lines (read group), preserved from input
	PGLines    []string          // @PG lines (program), preserved from input
}

// Reference represents a reference sequence (chromosome)
type Reference struct {
	ID   int32  // Reference ID
	Name string // Reference name
	Len  int32  // Reference length
}

// Record represents a single BAM alignment record
type Record struct {
	Name      string      // Query name
	Flags     uint16      // Bitwise flags
	RefID     int32       // Reference ID (-1 if unmapped)
	Pos       int32       // 0-based leftmost position
	MapQ      uint8       // Mapping quality
	Cigar     []CigarOp   // CIGAR operations
	MateRefID int32       // Mate reference ID
	MatePos   int32       // 0-based mate position
	TLen      int32       // Template length
	Seq       string      // Sequence
	Qual      []byte      // Quality scores (Phred+33)
	Aux       []*AuxField // Auxiliary fields
}

// CigarOp represents a single CIGAR operation
type CigarOp struct {
	Op  byte // Operation type (M, I, D, N, S, H, P, =, X)
	Len int  // Operation length
}

// CIGAR operation constants
const (
	CigarMatch     = 'M'
	CigarInsertion = 'I'
	CigarDeletion  = 'D'
	CigarSkip      = 'N' // Reference skip (intron)
	CigarSoftClip  = 'S'
	CigarHardClip  = 'H'
	CigarPadding   = 'P'
	CigarEqual     = '='
	CigarMismatch  = 'X'
)

// Flag constants
const (
	FlagPaired        = 1 << 0  // Template having multiple segments
	FlagProperPair    = 1 << 1  // Each segment properly aligned
	FlagUnmapped      = 1 << 2  // Segment unmapped
	FlagMateUnmapped  = 1 << 3  // Next segment unmapped
	FlagReverse       = 1 << 4  // SEQ being reverse complemented
	FlagMateReverse   = 1 << 5  // Next SEQ being reverse complemented
	FlagFirstInPair   = 1 << 6  // First segment
	FlagSecondInPair  = 1 << 7  // Last segment
	FlagSecondary     = 1 << 8  // Secondary alignment
	FlagQCFail        = 1 << 9  // Not passing quality controls
	FlagDuplicate     = 1 << 10 // PCR or optical duplicate
	FlagSupplementary = 1 << 11 // Supplementary alignment
)

// AuxField represents a BAM auxiliary field (tag-value pair)
type AuxField struct {
	Tag       string // 2-character tag
	Type      byte   // Value type (A, c, C, s, S, i, I, f, d, Z, H, B)
	ArrayType byte   // For array types, the element type (0 if not an array)
	Value     interface{}
}

// Aux type constants
const (
	AuxTypeChar   = 'A' // Printable character
	AuxTypeInt8   = 'c' // Signed 8-bit integer
	AuxTypeUInt8  = 'C' // Unsigned 8-bit integer
	AuxTypeInt16  = 's' // Signed 16-bit integer
	AuxTypeUInt16 = 'S' // Unsigned 16-bit integer
	AuxTypeInt32  = 'i' // Signed 32-bit integer
	AuxTypeUInt32 = 'I' // Unsigned 32-bit integer
	AuxTypeFloat  = 'f' // Single-precision floating point
	AuxTypeDouble = 'd' // Double-precision floating point
	AuxTypeString = 'Z' // Printable string
	AuxTypeHex    = 'H' // Byte array in Hex
	AuxTypeArray  = 'B' // Numeric array
)

// Reader reads BAM files
type Reader struct {
	header *Header
	r      *bgzip.Reader
}

// NewReader creates a new BAM reader
func NewReader(r io.Reader) (*Reader, error) {
	// Wrap in BGZF reader if not already
	bgzf, err := bgzip.NewReader(r)
	if err != nil {
		return nil, err
	}

	// Read and verify BAM magic
	magic := make([]byte, 4)
	_, err = bgzf.ReadFull(magic)
	if err != nil {
		return nil, fmt.Errorf("failed to read BAM magic: %w", err)
	}

	if string(magic) != "BAM\x01" {
		return nil, fmt.Errorf("invalid BAM magic: %v", magic)
	}

	// Read header
	header, err := readHeader(bgzf)
	if err != nil {
		return nil, err
	}

	return &Reader{
		header: header,
		r:      bgzf,
	}, nil
}

// Read reads a single BAM record
func (br *Reader) Read() (*Record, error) {
	return readRecord(br.r, br.header)
}

// VirtualOffset returns the BAI virtual offset of the next byte to be read
// from the underlying BGZF stream.
func (br *Reader) VirtualOffset() int64 {
	return br.r.VirtualOffset()
}

// Header returns the BAM header
func (br *Reader) Header() *Header {
	return br.header
}

// readHeader reads the BAM header
func readHeader(r *bgzip.Reader) (*Header, error) {
	// Read header text length (4 bytes)
	lenBuf := make([]byte, 4)
	_, err := r.ReadFull(lenBuf)
	if err != nil {
		return nil, fmt.Errorf("failed to read header text length: %w", err)
	}

	headerTextLen := int32(binary.LittleEndian.Uint32(lenBuf))
	if headerTextLen <= 0 {
		return nil, ErrInvalidHeader
	}

	// Read header text
	headerTextBytes := make([]byte, headerTextLen)
	_, err = r.ReadFull(headerTextBytes)
	if err != nil {
		return nil, fmt.Errorf("failed to read header text: %w", err)
	}

	// Verify we read the correct amount
	if len(headerTextBytes) != int(headerTextLen) {
		return nil, fmt.Errorf("header text length mismatch: expected %d, got %d", headerTextLen, len(headerTextBytes))
	}

	headerText := string(headerTextBytes)

	// Parse header
	header := &Header{
		OtherLines: make(map[string]string),
	}

	lines := strings.Split(headerText, "\n")
	for _, line := range lines {
		if len(line) == 0 {
			continue
		}
		fields := strings.Split(line, "\t")
		if len(fields) == 0 {
			continue
		}

		tag := fields[0]
		switch tag {
		case "@HD":
			for _, field := range fields[1:] {
				if strings.HasPrefix(field, "VN:") {
					header.Version = strings.TrimPrefix(field, "VN:")
				} else if strings.HasPrefix(field, "SO:") {
					header.SortOrder = strings.TrimPrefix(field, "SO:")
				}
			}
		case "@SQ":
			// Reference sequence
			ref := &Reference{}
			for _, field := range fields[1:] {
				if strings.HasPrefix(field, "SN:") {
					ref.Name = strings.TrimPrefix(field, "SN:")
				} else if strings.HasPrefix(field, "LN:") {
					fmt.Sscanf(strings.TrimPrefix(field, "LN:"), "%d", &ref.Len)
				}
			}
			if ref.Name != "" {
				ref.ID = int32(len(header.References))
				header.References = append(header.References, ref)
			}
		case "@RG":
			header.RGLines = append(header.RGLines, line)
		case "@PG":
			header.PGLines = append(header.PGLines, line)
		default:
			header.OtherLines[tag] = line
		}
	}

	// Read number of reference sequences
	nRefBuf := make([]byte, 4)
	_, err = r.ReadFull(nRefBuf)
	if err != nil {
		return nil, fmt.Errorf("failed to read reference count: %w", err)
	}

	nRef := int32(binary.LittleEndian.Uint32(nRefBuf))

	// Validate reference count matches parsed @SQ lines
	if nRef != int32(len(header.References)) {
		return nil, ErrInvalidHeader
	}

	// Read binary reference data (name_len + name + length for each reference)
	// We need to skip this even though we already parsed from text
	for i := int32(0); i < nRef; i++ {
		// Read name length
		nameLenBuf := make([]byte, 4)
		_, err := r.ReadFull(nameLenBuf)
		if err != nil {
			return nil, fmt.Errorf("failed to read ref name length: %w", err)
		}
		nameLen := int32(binary.LittleEndian.Uint32(nameLenBuf))

		// Read name
		nameBuf := make([]byte, nameLen)
		_, err = r.ReadFull(nameBuf)
		if err != nil {
			return nil, fmt.Errorf("failed to read ref name: %w", err)
		}

		// Read reference length
		lenBuf := make([]byte, 4)
		_, err = r.ReadFull(lenBuf)
		if err != nil {
			return nil, fmt.Errorf("failed to read ref length: %w", err)
		}
	}

	return header, nil
}

// readRecord reads a single BAM record
func readRecord(r *bgzip.Reader, header *Header) (*Record, error) {
	// Read block size
	blockSizeBuf := make([]byte, 4)
	_, err := r.ReadFull(blockSizeBuf)
	if err != nil {
		if err == io.EOF {
			return nil, io.EOF
		}
		return nil, fmt.Errorf("failed to read block size: %w", err)
	}

	blockSize := int32(binary.LittleEndian.Uint32(blockSizeBuf))

	if blockSize <= 0 {
		return nil, ErrInvalidRecord
	}

	// Read the rest of the record
	recordData := make([]byte, blockSize)
	_, err = r.ReadFull(recordData)
	if err != nil {
		return nil, fmt.Errorf("failed to read record data: %w", err)
	}

	// Parse record
	buf := recordData
	record := &Record{}

	// Read RefID
	if len(buf) < 4 {
		return nil, ErrInvalidRecord
	}
	record.RefID = int32(binary.LittleEndian.Uint32(buf[0:4]))
	buf = buf[4:]

	// Read Position
	if len(buf) < 4 {
		return nil, ErrInvalidRecord
	}
	record.Pos = int32(binary.LittleEndian.Uint32(buf[0:4]))
	buf = buf[4:]

	// Read Bin_mq_nl (bin: upper 16, mapq: next 8, l_qname: lower 8)
	if len(buf) < 4 {
		return nil, ErrInvalidRecord
	}
	binMQNL := binary.LittleEndian.Uint32(buf[0:4])
	buf = buf[4:]

	// Extract l_qname (lowest byte of bin_mq_nl)
	l_qname := int(binMQNL & 0xFF)
	record.MapQ = uint8((binMQNL >> 8) & 0xFF)

	// Read flag_nc (flag: upper 16, n_cigar_op: lower 16)
	if len(buf) < 4 {
		return nil, ErrInvalidRecord
	}
	flagNC := binary.LittleEndian.Uint32(buf[0:4])
	buf = buf[4:]

	// Extract flags and cigar count
	record.Flags = uint16(flagNC >> 16)
	nCigarOp := int(flagNC & 0xFFFF)

	// Read l_seq (4 bytes) from record core data.
	if len(buf) < 4 {
		return nil, ErrInvalidRecord
	}
	lSeq := int(binary.LittleEndian.Uint32(buf[0:4]))
	if lSeq < 0 {
		return nil, ErrInvalidRecord
	}
	buf = buf[4:]

	// Read Mate RefID
	if len(buf) < 4 {
		return nil, ErrInvalidRecord
	}
	record.MateRefID = int32(binary.LittleEndian.Uint32(buf[0:4]))
	buf = buf[4:]

	// Read Mate Position
	if len(buf) < 4 {
		return nil, ErrInvalidRecord
	}
	record.MatePos = int32(binary.LittleEndian.Uint32(buf[0:4]))
	buf = buf[4:]

	// Read Template Length
	if len(buf) < 4 {
		return nil, ErrInvalidRecord
	}
	record.TLen = int32(binary.LittleEndian.Uint32(buf[0:4]))
	buf = buf[4:]

	// Read Read Name (l_qname bytes, including null terminator)
	if len(buf) < l_qname {
		return nil, ErrInvalidRecord
	}
	record.Name = string(buf[0 : l_qname-1]) // Exclude null terminator
	buf = buf[l_qname:]

	// Read CIGAR operations
	record.Cigar = make([]CigarOp, nCigarOp)
	for i := 0; i < nCigarOp; i++ {
		if len(buf) < 4 {
			return nil, ErrInvalidRecord
		}
		cigarInt := binary.LittleEndian.Uint32(buf[0:4])
		buf = buf[4:]

		opLen := cigarInt >> 4
		opType := byte(cigarInt & 0xf)

		// Convert numeric op type to character
		opChar := cigarNumToChar(opType)
		record.Cigar[i] = CigarOp{Op: opChar, Len: int(opLen)}
	}

	// Read Sequence
	seqBytes := (lSeq + 1) / 2 // 2 bases per byte
	if len(buf) < seqBytes {
		return nil, ErrInvalidRecord
	}

	seqBuf := buf[0:seqBytes]
	buf = buf[seqBytes:]

	// Decode sequence (2-bit encoding)
	var seqBuilder strings.Builder
	seqBuilder.Grow(lSeq)
	for i := 0; i < lSeq; i++ {
		byteIdx := i / 2
		baseIdx := i % 2
		baseByte := seqBuf[byteIdx]

		var base byte
		if baseIdx == 0 {
			base = (baseByte >> 4) & 0xf
		} else {
			base = baseByte & 0xf
		}

		seqBuilder.WriteByte(seqBaseToChar(base))
	}
	record.Seq = seqBuilder.String()

	// Read Quality Scores
	if len(buf) < lSeq {
		return nil, ErrInvalidRecord
	}
	record.Qual = make([]byte, lSeq)
	copy(record.Qual, buf[0:lSeq])
	buf = buf[lSeq:]

	// Read Auxiliary Fields
	record.Aux = []*AuxField{}
	for len(buf) > 0 {
		if len(buf) < 3 {
			break
		}

		// Read Tag (2 bytes)
		tag := string(buf[0:2])
		buf = buf[2:]

		// Read Type (1 byte)
		auxType := buf[0]
		buf = buf[1:]

		aux := &AuxField{Tag: tag, Type: auxType}

		// Read Value based on type
		switch auxType {
		case AuxTypeChar:
			if len(buf) < 1 {
				return nil, ErrInvalidAux
			}
			aux.Value = string(buf[0:1])
			buf = buf[1:]
		case AuxTypeInt8:
			if len(buf) < 1 {
				return nil, ErrInvalidAux
			}
			aux.Value = int8(buf[0])
			buf = buf[1:]
		case AuxTypeUInt8:
			if len(buf) < 1 {
				return nil, ErrInvalidAux
			}
			aux.Value = buf[0]
			buf = buf[1:]
		case AuxTypeInt16:
			if len(buf) < 2 {
				return nil, ErrInvalidAux
			}
			aux.Value = int16(binary.LittleEndian.Uint16(buf[0:2]))
			buf = buf[2:]
		case AuxTypeUInt16:
			if len(buf) < 2 {
				return nil, ErrInvalidAux
			}
			aux.Value = binary.LittleEndian.Uint16(buf[0:2])
			buf = buf[2:]
		case AuxTypeInt32:
			if len(buf) < 4 {
				return nil, ErrInvalidAux
			}
			aux.Value = int32(binary.LittleEndian.Uint32(buf[0:4]))
			buf = buf[4:]
		case AuxTypeUInt32:
			if len(buf) < 4 {
				return nil, ErrInvalidAux
			}
			aux.Value = binary.LittleEndian.Uint32(buf[0:4])
			buf = buf[4:]
		case AuxTypeFloat:
			if len(buf) < 4 {
				return nil, ErrInvalidAux
			}
			bits := binary.LittleEndian.Uint32(buf[0:4])
			aux.Value = math.Float32frombits(bits)
			buf = buf[4:]
		case AuxTypeDouble:
			if len(buf) < 8 {
				return nil, ErrInvalidAux
			}
			bits := binary.LittleEndian.Uint64(buf[0:8])
			aux.Value = math.Float64frombits(bits)
			buf = buf[8:]
		case AuxTypeString:
			nullIdx := bytesIndexByte(buf, 0)
			if nullIdx < 0 {
				return nil, ErrInvalidAux
			}
			aux.Value = string(buf[0:nullIdx])
			buf = buf[nullIdx+1:]
		case AuxTypeHex:
			nullIdx := bytesIndexByte(buf, 0)
			if nullIdx < 0 {
				return nil, ErrInvalidAux
			}
			aux.Value = buf[0:nullIdx]
			buf = buf[nullIdx+1:]
		case AuxTypeArray:
			// Array type: type (1 byte), length (4 bytes), values
			if len(buf) < 5 {
				return nil, ErrInvalidAux
			}
			arrayType := buf[0]
			arrayLen := int(binary.LittleEndian.Uint32(buf[1:5]))
			buf = buf[5:]

			elemSize := auxArrayElemSize(arrayType)
			if elemSize == 0 {
				return nil, ErrInvalidAux
			}
			byteLen := arrayLen * elemSize
			if arrayLen < 0 || byteLen < 0 || len(buf) < byteLen {
				return nil, ErrInvalidAux
			}
			aux.Value = buf[0:byteLen]
			aux.ArrayType = arrayType // Store the array element type
			buf = buf[byteLen:]
		default:
			return nil, ErrInvalidAux
		}

		record.Aux = append(record.Aux, aux)
	}

	return record, nil
}

// GetAuxField returns the auxiliary field with the given tag, or nil if not found
func (r *Record) GetAuxField(tag string) *AuxField {
	for _, aux := range r.Aux {
		if aux.Tag == tag {
			return aux
		}
	}
	return nil
}

// IsPaired returns true if the record is paired
func (r *Record) IsPaired() bool {
	return r.Flags&FlagPaired != 0
}

// IsUnmapped returns true if the record is unmapped
func (r *Record) IsUnmapped() bool {
	return r.Flags&FlagUnmapped != 0
}

// IsMateUnmapped returns true if the mate is unmapped
func (r *Record) IsMateUnmapped() bool {
	return r.Flags&FlagMateUnmapped != 0
}

// IsReverse returns true if the record is reverse complemented
func (r *Record) IsReverse() bool {
	return r.Flags&FlagReverse != 0
}

// IsFirstInPair returns true if this is the first read in pair
func (r *Record) IsFirstInPair() bool {
	return r.Flags&FlagFirstInPair != 0
}

// IsSecondInPair returns true if this is the second read in pair
func (r *Record) IsSecondInPair() bool {
	return r.Flags&FlagSecondInPair != 0
}

// IsSecondary returns true if this is a secondary alignment
func (r *Record) IsSecondary() bool {
	return r.Flags&FlagSecondary != 0
}

// Helper functions

func bytesIndexByte(buf []byte, target byte) int {
	for i, b := range buf {
		if b == target {
			return i
		}
	}
	return -1
}

func auxArrayElemSize(t byte) int {
	switch t {
	case AuxTypeInt8, AuxTypeUInt8:
		return 1
	case AuxTypeInt16, AuxTypeUInt16:
		return 2
	case AuxTypeInt32, AuxTypeUInt32, AuxTypeFloat:
		return 4
	default:
		return 0
	}
}

func cigarNumToChar(n byte) byte {
	switch n {
	case 0:
		return CigarMatch
	case 1:
		return CigarInsertion
	case 2:
		return CigarDeletion
	case 3:
		return CigarSkip
	case 4:
		return CigarSoftClip
	case 5:
		return CigarHardClip
	case 6:
		return CigarPadding
	case 7:
		return CigarEqual
	case 8:
		return CigarMismatch
	default:
		return '?'
	}
}

func seqBaseToChar(b byte) byte {
	switch b {
	case 0:
		return '='
	case 1:
		return 'A'
	case 2:
		return 'C'
	case 4:
		return 'G'
	case 8:
		return 'T'
	case 15:
		return 'N'
	default:
		return 'N' // other IUPAC ambiguity codes
	}
}
