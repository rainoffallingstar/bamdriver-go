package bgzip

import (
	"bytes"
	"fmt"
	"io"
	"os"

	kpgzip "github.com/klauspost/compress/gzip"
)

// Writer writes BGZF compressed data
type Writer struct {
	w                 io.Writer
	buf               *bytes.Buffer
	blockSize         int
	closed            bool
	compressedWritten int64 // bytes written to the underlying file so far
}

// NewWriter creates a new BGZF writer
func NewWriter(path string) (*Writer, error) {
	f, err := os.Create(path)
	if err != nil {
		return nil, err
	}

	return &Writer{
		w:         f,
		buf:       bytes.NewBuffer(nil),
		blockSize: 0,
	}, nil
}

// Write writes data to the internal buffer, splitting across block boundaries
// so that each BGZF block contains at most BGZF_MAX_BLOCK_SIZE uncompressed bytes.
// htslib allocates exactly BGZF_MAX_BLOCK_SIZE bytes for decompression output;
// blocks that exceed this limit cause inflate() to return Z_BUF_ERROR.
func (z *Writer) Write(p []byte) (int, error) {
	if z.closed {
		return 0, fmt.Errorf("writer is closed")
	}

	total := 0
	for len(p) > 0 {
		// How many bytes we can still add to the current block
		space := BGZF_MAX_BLOCK_SIZE - z.blockSize
		chunk := p
		if len(chunk) > space {
			chunk = p[:space]
		}

		z.buf.Write(chunk)
		z.blockSize += len(chunk)
		p = p[len(chunk):]
		total += len(chunk)

		if z.blockSize >= BGZF_MAX_BLOCK_SIZE {
			if err := z.flushBlock(); err != nil {
				return total, err
			}
		}
	}

	return total, nil
}

// WriteByte writes a single byte
func (z *Writer) WriteByte(b byte) error {
	p := []byte{b}
	_, err := z.Write(p)
	return err
}

// flushBlock writes the current buffer as a BGZF block.
//
// It uses klauspost/compress/gzip so the DEFLATE payload is compatible with
// htslib's zlib-based BGZF reader.  After compression we patch the two-byte
// BSIZE field that lives inside the BC extra sub-field.
//
// Gzip block layout when FEXTRA is set with a 6-byte extra field:
//
//	Bytes  0-9   standard gzip header (ID1 ID2 CM FLG MTIME XFL OS)
//	Bytes 10-11  XLEN = 6  (little-endian)
//	Bytes 12-13  SI1='B' SI2='C'
//	Bytes 14-15  SLEN = 2  (little-endian)
//	Bytes 16-17  BSIZE placeholder  ← patched here
//	Bytes 18+    DEFLATE data, CRC32, ISIZE
func (z *Writer) flushBlock() error {
	if z.buf.Len() == 0 {
		return nil
	}

	uncompressed := z.buf.Bytes()

	// Compress using klauspost gzip with a BC extra-field placeholder.
	var compressed bytes.Buffer
	gw, err := kpgzip.NewWriterLevel(&compressed, kpgzip.DefaultCompression)
	if err != nil {
		return fmt.Errorf("failed to create gzip writer: %w", err)
	}
	// Pre-allocate BC sub-field: SI1 SI2 SLEN(2) BSIZE(2, will be patched)
	gw.Extra = []byte{
		BGZF_SI1, BGZF_SI2,
		byte(BGZF_SLEN), 0x00,
		0x00, 0x00, // BSIZE placeholder
	}

	if _, err := gw.Write(uncompressed); err != nil {
		return fmt.Errorf("failed to compress block: %w", err)
	}
	if err := gw.Close(); err != nil {
		return fmt.Errorf("failed to close gzip writer: %w", err)
	}

	block := compressed.Bytes()

	// Patch BSIZE = (total block size) - 1
	bsize := len(block) - 1
	block[16] = byte(bsize & 0xff)
	block[17] = byte((bsize >> 8) & 0xff)

	if _, err := z.w.Write(block); err != nil {
		return fmt.Errorf("failed to write BGZF block: %w", err)
	}
	z.compressedWritten += int64(len(block))

	z.buf.Reset()
	z.blockSize = 0

	return nil
}

// VirtualOffset returns the BAI virtual offset of the next byte to be written.
// compressedWritten is the block start in the file; buf.Len() is the offset within
// the current (not-yet-flushed) uncompressed block.
func (z *Writer) VirtualOffset() int64 {
	return (z.compressedWritten << 16) | int64(z.buf.Len())
}

// Close flushes remaining data, writes the BGZF EOF sentinel block, and
// closes the underlying file.
func (z *Writer) Close() error {
	if z.closed {
		return nil
	}

	// Flush any remaining data
	if z.buf.Len() > 0 {
		if err := z.flushBlock(); err != nil {
			return err
		}
	}

	// Write BGZF EOF sentinel block (SAM spec §4.1).
	// This is a fixed 28-byte empty BGZF block required by all BGZF readers.
	eofBlock := []byte{
		0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00,
		0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
		0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x00, 0x00, 0x00, 0x00,
	}
	if _, err := z.w.Write(eofBlock); err != nil {
		return fmt.Errorf("failed to write BGZF EOF block: %w", err)
	}

	z.closed = true

	// Close underlying writer if it implements io.Closer
	if closer, ok := z.w.(io.Closer); ok {
		return closer.Close()
	}

	return nil
}
