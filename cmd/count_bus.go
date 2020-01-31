package cmd

import (
	"bufio"
	"encoding/binary"
	"fmt"
	"io"
	"log"
	"os"
	"strconv"
	"strings"

	"github.com/spf13/cobra"
)

var busFile string
var geneMapFile string
var ecMapFile string
var txNamesFile string
var cellsFile string

// var outputFile string
var outputDir string
var includeMultigeneECs bool
var calcReadCount bool

func init() {
	// Sorted BUS file with corrected barcodes (output.correct.sort.bus)
	// TODO: File with transcript compatibility counts (matrix.tcc.mtx)
	countCmd.Flags().StringVarP(&busFile, "bus", "b", "", "Input BUS file")
	// File with gene names (transcripts_to_genes.txt)
	countCmd.Flags().StringVarP(&geneMapFile, "genemap", "g", "", "File for mapping transcripts to genes")
	// File with equivalence classes (matrix.ec)
	countCmd.Flags().StringVarP(&ecMapFile, "ecmap", "e", "", "File for mapping equivalence classes to transcripts")
	// File with transcripts names (transcripts.txt)
	countCmd.Flags().StringVarP(&txNamesFile, "txnames", "t", "", "File with names of transcripts")
	// File with barcodes (matrix.cells)
	countCmd.Flags().StringVarP(&cellsFile, "cells", "c", "", "List of barcodes")

	// Define output
	// Output file
	// countCmd.Flags().StringVarP(&outputFile, "output", "o", "genes_bus_matrix.mtx", "File to write output to")
	// Output directory
	countCmd.Flags().StringVarP(&outputDir, "output", "o", "echinus_out/", "Directory to write output files to")

	// Define boolean flags
	// Note: currently only unique mappers are used to calculate gene abundances
	countCmd.Flags().BoolVarP(&includeMultigeneECs, "multimapping", "m", false, "Include bus records that pseudoalign to multiple genes")
	// Count reads rather than molecules
	countCmd.Flags().BoolVarP(&calcReadCount, "reads", "r", false, "Count reads instead of unique molecules")

	rootCmd.AddCommand(countCmd)
}

var countCmd = &cobra.Command{
	Use:   "count",
	Short: "Get the plain text version of a busfile",
	Long:  "Convert a busfile to a plain text barcode-UMI-EC-count file",
	Run: func(cmd *cobra.Command, args []string) {
		// File with transcript compatibility counts (matrix.tcc.mtx)
		Count(busFile)
	},
}

func sliceAtoi(xs []string) (xi []int, err error) {
	xi = make([]int, 0, len(xs))
	for _, s := range xs {
		i, err := strconv.Atoi(s)
		if err != nil {
			return xi, err
		}
		xi = append(xi, i)
	}
	return
}

// Read a plain text file with one column
func readPlainFile(fileName string) (fileContent []string, err error) {
	file, err := os.Open(fileName)
	if err != nil {
		return
	}
	defer file.Close()

	fileScanner := bufio.NewScanner(file)
	for fileScanner.Scan() {
		fileContent = append(fileContent, fileScanner.Text())
	}
	return
}

// Read a plain text file with multiple columns, first column being unique
func readMapFile(fileName string) (fileMap map[string][]string, err error) {
	fileMap = make(map[string][]string)

	file, err := os.Open(fileName)
	if err != nil {
		return
	}
	defer file.Close()

	fileScanner := bufio.NewScanner(file)
	var key string
	var line, values []string
	for fileScanner.Scan() {
		line = strings.Split(fileScanner.Text(), "\t")

		key = line[0]
		values = line[1:len(line)]
		fileMap[key] = values
	}
	return
}

// Return a unique subset of the string slice provided.
func unique(xs []string) (uniq []string) {
	uniq = make([]string, 0, len(xs))
	seen := make(map[string]bool)

	for _, val := range xs {
		if _, ok := seen[val]; !ok {
			seen[val] = true
			uniq = append(uniq, val)
		}
	}

	return
}

func convertTrxToGenes(ts []int, trx []string, geneMap map[string][]string) (gs []string) {
	for _, ti := range ts {
		gs = append(gs, geneMap[trx[ti]][1])
	}
	gs = unique(gs)
	return
}

func index(xs []string, x string) int {
	for k, v := range xs {
		if v == x {
			return k
		}
	}
	return -1
}

// Count molecules or reads for every cell and every gene in the BUS file
func Count(busFile string) {

	// Read input files
	ec, err := os.Open(ecMapFile)
	if err != nil {
		log.Fatal(err)
	}
	defer ec.Close()

	txNames, err := readPlainFile(txNamesFile)
	if err != nil {
		log.Fatal(err)
	}

	geneMap, err := readMapFile(geneMapFile)
	if err != nil {
		log.Fatal(err)
	}
	var genes []string
	for _, v := range geneMap {
		genes = append(genes, v[1])
	}
	genes = unique(genes)

	cells, err := readPlainFile(cellsFile)
	if err != nil {
		log.Fatal(err)
	}

	ecScanner := bufio.NewScanner(ec)
	ecGeneMap := make(map[int32][]string)

	// Define structure for counts
	cellGeneCounts := make([][]uint32, len(cells))
	for cell := range cellGeneCounts {
		cellGeneCounts[cell] = make([]uint32, len(genes))
	}

	for ecScanner.Scan() {
		ecAndTranscripts := strings.Split(ecScanner.Text(), "\t")

		ecIndexInt, err := strconv.Atoi(ecAndTranscripts[0])
		if err != nil {
			log.Fatal(err)
		}
		ecIndex := int32(ecIndexInt)

		ecTrx, err := sliceAtoi(strings.Split(ecAndTranscripts[1], ","))
		if err != nil {
			log.Fatal(err)
		}

		// Create EC to gene mapping
		if includeMultigeneECs {
			ecGenes := convertTrxToGenes(ecTrx, txNames, geneMap)
			ecGeneMap[ecIndex] = ecGenes
		} else {
			if len(ecTrx) == 1 {
				ecGenes := convertTrxToGenes(ecTrx, txNames, geneMap)
				ecGeneMap[ecIndex] = ecGenes
			}
		}
	}

	// Create gene to EC mapping
	// geneEcMap := reverseMap(ecGeneMap)

	// Parse counts matrix
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
	fmt.Printf("%+v\n", fileHeader)

	// Read BUS data
	var data BUSData
	var lineIndex int
	var nCounts uint64

	// Count TCCs for every gene in every cell
	for {
		err := binary.Read(busReader, binary.LittleEndian, &data)
		if err != nil {
			if err == io.EOF {
				break
			}
			log.Fatal(err)
		}

		if gene, ok := ecGeneMap[data.EC]; ok {
			cell := binaryToString(data.Barcode, fileHeader.BCLen)
			cellInt := index(cells, cell)

			// Note: this only works when a single gene corresponds to a single EC
			// TODO: implement counting with multi-mapping reads
			geneInt := index(genes, gene[0])

			// Consider an option to add up counts
			if calcReadCount {
				cellGeneCounts[cellInt][geneInt] += data.Count
				nCounts += uint64(data.Count)
			} else {
				cellGeneCounts[cellInt][geneInt]++
				nCounts++
			}
		}
		lineIndex++
		fmt.Printf("\rProcessed line %d, recorded %d counts", lineIndex, nCounts)
	}

	// Save the matrix file
	// outMtx, err := os.Create("genes_bus.mtx")
	outMtx, err := os.Create(outputDir + "/matrix.mtx")
	if err != nil {
		log.Fatal(err)
	}
	defer outMtx.Close()

	// Get total number of elements
	var nElements int // Store total number of elements
	for _, cellCounts := range cellGeneCounts {
		for _, count := range cellCounts {
			if count > 0 {
				nElements++
			}
		}
	}

	// Matrix header
	mtx := bufio.NewWriter(outMtx)
	_, err = mtx.WriteString("%%MatrixMarket matrix coordinate real general\n")
	if err != nil {
		log.Fatal(err)
	}

	_, err = mtx.WriteString(strconv.Itoa(len(cells)) + "\t" + strconv.Itoa(len(genes)) + "\t" + strconv.Itoa(nElements) + "\n")
	if err != nil {
		log.Fatal(err)
	}

	var count int
	for cellIndex := range cells {
		for geneIndex := range genes {
			count = int(cellGeneCounts[cellIndex][geneIndex])
			if count > 0 {
				_, err = mtx.WriteString(strconv.Itoa(cellIndex+1) + "\t" + strconv.Itoa(geneIndex+1) + "\t" + strconv.Itoa(count) + "\n")
				if err != nil {
					log.Fatal(err)
				}
			}
		}
	}
	mtx.Flush()

	// Write barcodes.tsv
	outBarcodesTSV, err := os.Create(outputDir + "/barcodes.tsv")
	if err != nil {
		log.Fatal(err)
	}
	defer outBarcodesTSV.Close()
	barcodesTSV := bufio.NewWriter(outBarcodesTSV)

	for _, cellBarcode := range cells {
		_, err = barcodesTSV.WriteString(cellBarcode + "\n")
		if err != nil {
			log.Fatal(err)
		}
	}
	barcodesTSV.Flush()

	// Write genes.tsv
	outGenesTSV, err := os.Create(outputDir + "/genes.tsv")
	if err != nil {
		log.Fatal(err)
	}
	defer outGenesTSV.Close()
	genesTSV := bufio.NewWriter(outGenesTSV)

	for _, geneName := range genes {
		_, err = genesTSV.WriteString(geneName + "\n")
		if err != nil {
			log.Fatal(err)
		}
	}
	genesTSV.Flush()
}
