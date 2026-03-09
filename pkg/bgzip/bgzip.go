// Package bgzip implements BGZF (Blocked GNU Zip Format) compression and decompression.
// BGZF is a variant of gzip that uses blocks of up to 64KB uncompressed data,
// allowing for random access to compressed BAM files.
package bgzip

import (
	"bytes"
	"errors"
	"fmt"
	"io"
	"os"

	kpgzip "github.com/klauspost/compress/gzip"
)

var (
	ErrHeader      = errors.New("invalid BGZF header")
	ErrTruncated   = errors.New("truncated BGZF block")
	ErrExtraField  = errors.New("invalid BGZF extra field")
	ErrBlockSize   = errors.New("invalid block size")
	ErrNoBGZF      = errors.New("not a BGZF file")
	ErrBadChecksum = errors.New("BGZF checksum mismatch")
)

// BGZF magic numbers and constants
const (
	BGZF_ID1            = 31        // First byte of gzip magic
	BGZF_ID2            = 139       // Second byte of gzip magic
	BGZF_CM_DEFLATE     = 8         // Compression method
	BGZF_FLG_FEXTRA     = 4         // Extra field flag
	BGZF_XLEN           = 6         // Length of extra field
	BGZF_SI1            = 66        // 'B'
	BGZF_SI2            = 67        // 'C'
	BGZF_SLEN           = 2         // Subfield length
	BGZF_MAX_BLOCK_SIZE = 64 * 1024 // Maximum uncompressed block size
)

// Block represents a single BGZF block
type Block struct {
	Data       []byte // Decompressed data
	Offset     int64  // File offset of this block
	UncompSize int    // Uncompressed size
}

// Reader reads BGZF compressed data
type Reader struct {
	r   io.Reader
	buf []byte // Buffer for decompressed data
	pos int    // Current position in buffer
	err error
	eof bool

	// Track compressed file offset for indexing
	compressedOffset int64 // Current compressed (file) offset
	lastBlockSize    int   // Last block size read
}

// NewReader creates a new BGZF reader
func NewReader(r io.Reader) (*Reader, error) {
	bg := &Reader{
		r:   r,
		buf: nil,
		pos: 0,
	}
	return bg, nil
}

// Read reads decompressed data from the BGZF stream
func (bg *Reader) Read(p []byte) (int, error) {
	if bg.err != nil {
		return 0, bg.err
	}

	// If buffer is empty or exhausted, read next block
	if bg.buf == nil || bg.pos >= len(bg.buf) {
		block, err := bg.readBlock()
		if err != nil {
			if err == io.EOF {
				// This is the real end of file
				bg.eof = true
				return 0, io.EOF
			}
			bg.err = err
			return 0, err
		}

		bg.buf = block.Data
		bg.pos = 0

		// BGZF EOF marker is a valid empty block (0 decompressed bytes).
		// Return io.EOF now so callers don't loop on (0, nil).
		if len(bg.buf) == 0 {
			bg.eof = true
			return 0, io.EOF
		}
	}

	// Copy from buffer to p
	n := copy(p, bg.buf[bg.pos:])
	bg.pos += n

	return n, nil
}

// readBlock reads and decompresses a single BGZF block
func (bg *Reader) readBlock() (*Block, error) {
	// Save the compressed offset of this block's start before reading any bytes.
	blockStart := bg.compressedOffset

	// Read standard gzip header (10 bytes)
	stdHeader := make([]byte, 10)
	n, err := io.ReadFull(bg.r, stdHeader)
	if err != nil {
		if err == io.EOF {
			return nil, io.EOF
		}
		return nil, fmt.Errorf("failed to read standard header: %w", err)
	}
	if n < 10 {
		return nil, ErrTruncated
	}

	// Validate gzip magic
	if stdHeader[0] != BGZF_ID1 || stdHeader[1] != BGZF_ID2 {
		return nil, ErrNoBGZF
	}

	// Validate compression method
	if stdHeader[2] != BGZF_CM_DEFLATE {
		return nil, ErrHeader
	}

	// Check for extra field flag
	if stdHeader[3] != BGZF_FLG_FEXTRA {
		return nil, ErrNoBGZF // Not BGZF
	}

	// Read XLEN (extra field length, 2 bytes, little-endian)
	xlenBytes := make([]byte, 2)
	n, err = io.ReadFull(bg.r, xlenBytes)
	if err != nil {
		return nil, fmt.Errorf("failed to read XLEN: %w", err)
	}
	xlen := int(xlenBytes[0]) | int(xlenBytes[1])<<8

	// Read the extra field
	extraField := make([]byte, xlen)
	n, err = io.ReadFull(bg.r, extraField)
	if err != nil {
		return nil, fmt.Errorf("failed to read extra field: %w", err)
	}
	if n < xlen {
		return nil, ErrTruncated
	}

	// Parse BGZF BC subfield from extra field
	foundBC := false
	var bsize int
	for i := 0; i <= xlen-6; i++ {
		if extraField[i] == BGZF_SI1 && extraField[i+1] == BGZF_SI2 {
			slen := int(extraField[i+2]) | int(extraField[i+3])<<8
			if slen != BGZF_SLEN {
				continue
			}
			bsize = int(extraField[i+4]) | int(extraField[i+5])<<8
			foundBC = true
			break
		}
	}

	if !foundBC {
		return nil, ErrNoBGZF
	}

	// Calculate total BGZF block size
	// BSIZE + 1 = total block size
	totalBlockSize := bsize + 1

	// Read the entire rest of the BGZF block (compressed data + 8-byte trailer)
	// We've already read: 10 (std header) + 2 (XLEN) + xlen (extra field) = 12 + xlen bytes
	// Remaining: totalBlockSize - (12 + xlen)
	remainingInBlock := totalBlockSize - (12 + xlen)
	restOfBlock := make([]byte, remainingInBlock)
	n, err = io.ReadFull(bg.r, restOfBlock)
	if err != nil {
		return nil, fmt.Errorf("failed to read rest of BGZF block (need %d bytes): %w", remainingInBlock, err)
	}
	if n < remainingInBlock {
		return nil, ErrTruncated
	}

	// Update compressed offset tracking now that all bytes of this block are read.
	bg.compressedOffset += int64(totalBlockSize)
	bg.lastBlockSize = totalBlockSize

	// Reconstruct complete gzip stream
	gzipStream := bytes.NewBuffer(nil)
	gzipStream.Write(stdHeader)   // 10 bytes
	gzipStream.Write(xlenBytes)   // 2 bytes
	gzipStream.Write(extraField)  // xlen bytes
	gzipStream.Write(restOfBlock) // compressed data + 8-byte trailer

	// Use klauspost gzip reader for decompression
	gzipReader, err := kpgzip.NewReader(gzipStream)
	if err != nil {
		return nil, fmt.Errorf("failed to create gzip reader: %w", err)
	}
	defer gzipReader.Close()

	// Read decompressed data
	decompressed, err := io.ReadAll(gzipReader)
	if err != nil {
		return nil, fmt.Errorf("decompression failed: %w", err)
	}

	return &Block{
		Data:       decompressed,
		Offset:     blockStart,
		UncompSize: len(decompressed),
	}, nil
}

// VirtualOffset returns the BAI virtual offset of the next byte to be read.
// Call before Read()/ReadFull() to get voffBeg, after to get voffEnd.
// voff = (compressed_block_start << 16) | offset_within_uncompressed_block
func (bg *Reader) VirtualOffset() int64 {
	if bg.buf == nil || bg.pos >= len(bg.buf) {
		// Buffer exhausted; next byte is at the start of the next block.
		return bg.compressedOffset << 16
	}
	blockStart := bg.compressedOffset - int64(bg.lastBlockSize)
	return (blockStart << 16) | int64(bg.pos)
}

// WriteTo writes the decompressed data to w
func (bg *Reader) WriteTo(w io.Writer) (int64, error) {
	// Don't use io.Copy as it causes infinite recursion
	// Instead, implement our own copy loop
	var written int64
	buf := make([]byte, 32*1024) // 32KB buffer

	for {
		n, err := bg.Read(buf)
		if n > 0 {
			nw, errw := w.Write(buf[:n])
			written += int64(nw)
			if errw != nil {
				return written, errw
			}
			if nw != n {
				return written, io.ErrShortWrite
			}
		}
		if err != nil {
			if err == io.EOF {
				return written, nil
			}
			return written, err
		}
	}
}

// ReadAll reads and decompresses all data from r
func ReadAll(r io.Reader) ([]byte, error) {
	bg, err := NewReader(r)
	if err != nil {
		return nil, err
	}

	var buf bytes.Buffer
	_, err = io.Copy(&buf, bg)
	if err != nil {
		return nil, err
	}

	return buf.Bytes(), nil
}

// ReadFull is like io.ReadFull but handles BGZF block boundaries correctly.
// It ensures we read exactly len(p) bytes, handling the case where a single
// Read() call may return fewer bytes due to BGZF block boundaries.
func (bg *Reader) ReadFull(p []byte) (int, error) {
	totalRead := 0
	for totalRead < len(p) {
		n, err := bg.Read(p[totalRead:])
		if err != nil {
			if err == io.EOF && totalRead > 0 {
				return totalRead, io.ErrUnexpectedEOF
			}
			return totalRead, err
		}
		if n == 0 {
			return totalRead, io.ErrNoProgress
		}
		totalRead += n
	}
	return totalRead, nil
}

// ReadFile reads and decompresses a BGZF file
func ReadFile(path string) ([]byte, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	return ReadAll(f)
}
