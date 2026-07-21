package bgzip

import (
	"bytes"
	"io"
	"testing"
)

func FuzzReaderMalformedInput(fuzz *testing.F) {
	fuzz.Add([]byte{})
	fuzz.Add([]byte{0x1f, 0x8b})
	fuzz.Add([]byte{
		31, 139, 8, 4, 0, 0, 0, 0, 0, 255,
		6, 0, 'B', 'C', 2, 0, 27, 0,
		3, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	})

	fuzz.Fuzz(func(t *testing.T, compressedData []byte) {
		reader, err := NewReader(bytes.NewReader(compressedData))
		if err != nil {
			return
		}
		buffer := make([]byte, 4096)
		for totalRead := 0; totalRead <= BGZF_MAX_BLOCK_SIZE*2; {
			bytesRead, readErr := reader.Read(buffer)
			totalRead += bytesRead
			if readErr != nil {
				if readErr == io.EOF {
					return
				}
				return
			}
			if bytesRead == 0 {
				t.Fatal("Reader returned no progress without an error")
			}
		}
	})
}
