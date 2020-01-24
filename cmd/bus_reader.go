package cmd

import (
	"bufio"
	"encoding/binary"
	"fmt"
	"io"
	"log"
	"os"
)

const busFormatVersion = 1

// BUSTranscript is to store a single transcript name and length
type BUSTranscript struct {
	Name             string
	TranscriptLength uint32
}

// BUSHeader is a header of a BUS file
// Transcripts and ECs are populated from external files
type BUSHeader struct {
	Text        string
	Version     uint32
	BCLen       uint32
	UMILen      uint32
	Transcripts []BUSTranscript
	ECs         [][]int32
}

// BUSData is a single record in a BUS file
type BUSData struct {
	Barcode uint64
	UMI     uint64
	EC      int32
	Count   uint32
	Flags   uint32
	Pad     uint32
}

func parseHeader(reader *bufio.Reader) (header BUSHeader) {
	header = BUSHeader{}

	buf := make([]byte, 4)

	reader.Read(buf)

	magic := string(buf)
	if magic != "BUS\000" {
		log.Fatal("File does not begin with the magic string BUS\000\n")
	}

	reader.Read(buf)
	header.Version = binary.LittleEndian.Uint32(buf)

	reader.Read(buf)
	header.BCLen = binary.LittleEndian.Uint32(buf)

	reader.Read(buf)
	header.UMILen = binary.LittleEndian.Uint32(buf)

	reader.Read(buf)
	tlen := binary.LittleEndian.Uint32(buf)

	tbuf := make([]byte, tlen)
	_, _ = reader.Read(tbuf)
	header.Text = string(tbuf)

	return
}

func binaryToString(x uint64, n uint32) string {
	var c rune
	cs := make([]rune, n)
	sh := uint(n - 1)
	var i uint32
	for i = 0; i < n; i++ {
		// c = 'N'
		switch (x >> (2 * sh)) & 0x03 {
		case 0x00:
			c = 'A'
		case 0x01:
			c = 'C'
		case 0x02:
			c = 'G'
		case 0x03:
			c = 'T'
		default:
			c = 'N'
		}
		sh--
		cs[i] = c
	}
	return string(cs)
}

// ReadBUS reads the BUS file and output its text representation
func ReadBUS(busFile string) error {

	bus, err := os.Open(busFile)
	if err != nil {
		log.Fatal(err)
	}
	defer bus.Close()

	busReader := bufio.NewReader(bus)

	fileHeader := parseHeader(busReader)
	if fileHeader.Version != busFormatVersion {
		log.Fatal("File version does not match the accepted BUS version")
	}

	// Read BUS data
	var data BUSData

	for {
		err := binary.Read(busReader, binary.LittleEndian, &data)
		if err != nil {
			if err == io.EOF {
				break
			}
			log.Fatal(err)
		}
		fmt.Println(binaryToString(data.Barcode, fileHeader.BCLen) + "\t" + binaryToString(data.UMI, fileHeader.UMILen) + "\t" + fmt.Sprint(data.EC) + "\t" + fmt.Sprint(data.Count))
	}

	fmt.Printf("%+v\n", fileHeader)

	return nil
}
