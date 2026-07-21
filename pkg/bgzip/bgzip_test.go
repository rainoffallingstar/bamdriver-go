package bgzip

import (
	"bytes"
	"encoding/binary"
	"errors"
	"io"
	"testing"
)

func TestWriterRoundTripIncompressibleData(t *testing.T) {
	input := makeDeterministicData(maximumUncompressedPayload * 3)
	var compressed bytes.Buffer
	writer := &Writer{w: &compressed, buf: bytes.NewBuffer(nil)}

	written, err := writer.Write(input)
	if err != nil {
		t.Fatalf("Write: %v", err)
	}
	if written != len(input) {
		t.Fatalf("written = %d, want %d", written, len(input))
	}
	if err := writer.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	assertValidBlockSizes(t, compressed.Bytes())
	output, err := ReadAll(bytes.NewReader(compressed.Bytes()))
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if !bytes.Equal(output, input) {
		t.Fatal("round-trip data mismatch")
	}
}

func TestReaderRejectsPartialHeader(t *testing.T) {
	reader, err := NewReader(bytes.NewReader([]byte{BGZF_ID1, BGZF_ID2, BGZF_CM_DEFLATE}))
	if err != nil {
		t.Fatalf("NewReader: %v", err)
	}

	buffer := make([]byte, 1)
	_, err = reader.Read(buffer)
	if !errors.Is(err, ErrTruncated) {
		t.Fatalf("Read error = %v, want ErrTruncated", err)
	}
}

func TestReaderRejectsDeclaredBlockSmallerThanHeader(t *testing.T) {
	block := []byte{
		BGZF_ID1, BGZF_ID2, BGZF_CM_DEFLATE, BGZF_FLG_FEXTRA,
		0, 0, 0, 0, 0, 255,
		6, 0,
		BGZF_SI1, BGZF_SI2, BGZF_SLEN, 0,
		10, 0,
	}
	reader, err := NewReader(bytes.NewReader(block))
	if err != nil {
		t.Fatalf("NewReader: %v", err)
	}

	_, err = reader.Read(make([]byte, 1))
	if !errors.Is(err, ErrBlockSize) {
		t.Fatalf("Read error = %v, want ErrBlockSize", err)
	}
}

func TestReaderRejectsMalformedExtraSubfield(t *testing.T) {
	block := []byte{
		BGZF_ID1, BGZF_ID2, BGZF_CM_DEFLATE, BGZF_FLG_FEXTRA,
		0, 0, 0, 0, 0, 255,
		6, 0,
		'X', 'Y', 10, 0, 0, 0,
	}
	reader, err := NewReader(bytes.NewReader(block))
	if err != nil {
		t.Fatalf("NewReader: %v", err)
	}

	_, err = reader.Read(make([]byte, 1))
	if !errors.Is(err, ErrExtraField) {
		t.Fatalf("Read error = %v, want ErrExtraField", err)
	}
}

func TestWriteFullHandlesShortWrites(t *testing.T) {
	writer := &limitedWriter{maximumPerWrite: 3}
	input := []byte("abcdefghij")
	if err := writeFull(writer, input); err != nil {
		t.Fatalf("writeFull: %v", err)
	}
	if !bytes.Equal(writer.buffer.Bytes(), input) {
		t.Fatalf("written data = %q, want %q", writer.buffer.Bytes(), input)
	}
}

func TestWriteFullRejectsNoProgress(t *testing.T) {
	err := writeFull(zeroWriter{}, []byte("data"))
	if !errors.Is(err, io.ErrShortWrite) {
		t.Fatalf("writeFull error = %v, want io.ErrShortWrite", err)
	}
}

func assertValidBlockSizes(t *testing.T, compressed []byte) {
	t.Helper()
	for offset := 0; offset < len(compressed); {
		if len(compressed)-offset < 18 {
			t.Fatalf("truncated block header at offset %d", offset)
		}
		blockSize := int(binary.LittleEndian.Uint16(compressed[offset+16:offset+18])) + 1
		if blockSize > BGZF_MAX_BLOCK_SIZE {
			t.Fatalf("block size = %d, exceeds %d", blockSize, BGZF_MAX_BLOCK_SIZE)
		}
		if blockSize < 18 || offset+blockSize > len(compressed) {
			t.Fatalf("invalid block size %d at offset %d", blockSize, offset)
		}
		offset += blockSize
	}
}

func makeDeterministicData(length int) []byte {
	data := make([]byte, length)
	state := uint32(0x9e3779b9)
	for index := range data {
		state ^= state << 13
		state ^= state >> 17
		state ^= state << 5
		data[index] = byte(state)
	}
	return data
}

type limitedWriter struct {
	buffer          bytes.Buffer
	maximumPerWrite int
}

func (writer *limitedWriter) Write(data []byte) (int, error) {
	if len(data) > writer.maximumPerWrite {
		data = data[:writer.maximumPerWrite]
	}
	return writer.buffer.Write(data)
}

type zeroWriter struct{}

func (zeroWriter) Write([]byte) (int, error) {
	return 0, nil
}
