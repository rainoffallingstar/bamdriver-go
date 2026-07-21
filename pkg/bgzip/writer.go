package bgzip

import (
	"bytes"
	"errors"
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

const maximumUncompressedPayload = 65280

// Write writes data to the internal buffer, splitting across block boundaries.
// The compressed BGZF member, including headers and trailer, must fit in a
// 16-bit BSIZE field. Using a conservative payload leaves enough room for
// DEFLATE overhead even when the input is incompressible.
func (z *Writer) Write(p []byte) (int, error) {
	if z.closed {
		return 0, fmt.Errorf("writer is closed")
	}

	total := 0
	for len(p) > 0 {
		space := maximumUncompressedPayload - z.blockSize
		chunk := p
		if len(chunk) > space {
			chunk = p[:space]
		}

		z.buf.Write(chunk)
		z.blockSize += len(chunk)
		p = p[len(chunk):]
		total += len(chunk)

		if z.blockSize >= maximumUncompressedPayload {
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
	if len(block) > BGZF_MAX_BLOCK_SIZE {
		return ErrBlockSize
	}
	if len(block) < 18 {
		return ErrBlockSize
	}

	// Patch BSIZE = (total block size) - 1.
	bsize := len(block) - 1
	block[16] = byte(bsize & 0xff)
	block[17] = byte((bsize >> 8) & 0xff)

	if err := writeFull(z.w, block); err != nil {
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
	z.closed = true

	var closeErrors []error
	if z.buf.Len() > 0 {
		if err := z.flushBlock(); err != nil {
			closeErrors = append(closeErrors, err)
		}
	}

	if len(closeErrors) == 0 {
		eofBlock := []byte{
			0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00,
			0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
			0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
			0x00, 0x00, 0x00, 0x00,
		}
		if err := writeFull(z.w, eofBlock); err != nil {
			closeErrors = append(closeErrors, fmt.Errorf("failed to write BGZF EOF block: %w", err))
		}
	}

	if closer, ok := z.w.(io.Closer); ok {
		if err := closer.Close(); err != nil {
			closeErrors = append(closeErrors, err)
		}
	}

	return errors.Join(closeErrors...)
}

func writeFull(writer io.Writer, data []byte) error {
	for len(data) > 0 {
		written, err := writer.Write(data)
		if err != nil {
			return err
		}
		if written <= 0 {
			return io.ErrShortWrite
		}
		data = data[written:]
	}
	return nil
}
