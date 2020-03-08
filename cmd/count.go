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
var tccFile string
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
	countCmd.Flags().StringVarP(&busFile, "bus", "b", "", "Input BUS file")
	// Alternative: file with transcript compatibility counts (matrix.tcc.mtx)
	countCmd.Flags().StringVarP(&tccFile, "tcc", "x", "", "Input .tcc.mtx matrix file")
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
	Short: "Generate count matrix per gene",
	Long:  "Count molecules (or reads) per gene in a BUS file (or a TCC matrix file)",
	Run: func(cmd *cobra.Command, args []string) {
		// If a file with transcript compatibility counts (matrix.tcc.mtx) is provided
		if tccFile != "" {
			CountMtx(tccFile)
		} else {
			// Use CountBus by default
			CountBus(busFile)
		}
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

// Return a unique subset of the string slice provided, with unique elements indices.
func uniqueWithIndices(xs []string) (uniq []string, uniqIndices []int) {
	uniq = make([]string, 0, len(xs))
	uniqIndices = make([]int, 0, len(xs))
	seen := make(map[string]bool)

	for i, val := range xs {
		if _, ok := seen[val]; !ok {
			seen[val] = true
			uniq = append(uniq, val)
			uniqIndices = append(uniqIndices, i)
		}
	}

	return
}

// Return a subset of the string slice provided using provided indices.
func subset(xs []string, uniqIndices []int) (uniq []string) {
	uniq = make([]string, 0, len(xs))

	for _, val := range uniqIndices {
		uniq = append(uniq, xs[val])
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

func convertTrxToGenes2col(ts []int, trx []string, geneMap map[string][]string) (gs []string) {
	for _, ti := range ts {
		gs = append(gs, geneMap[trx[ti]][0])
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

// CountBus count molecules or reads for every cell and every gene in the BUS file
func CountBus(busFile string) {

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
	var geneIDs []string
	var geneIndices []int
	var genes []string
	for _, v := range geneMap {
		geneIDs = append(geneIDs, v[0])
		genes = append(genes, v[1])
	}
	geneIDs, geneIndices = uniqueWithIndices(geneIDs)
	genes = subset(genes, geneIndices)

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
	fmt.Printf("\n")

	// Create output dir if does not exist
	if _, err := os.Stat(outputDir); os.IsNotExist(err) {
		os.Mkdir(outputDir, 0700)
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

	_, err = mtx.WriteString(strconv.Itoa(len(genes)) + "\t" + strconv.Itoa(len(cells)) + "\t" + strconv.Itoa(nElements) + "\n")
	if err != nil {
		log.Fatal(err)
	}

	var count int
	var cellCounts []uint32
	for cellIndex := range cells {
		cellCounts = cellGeneCounts[cellIndex]
		for geneIndex := range genes {
			count = int(cellCounts[geneIndex])
			if count > 0 {
				_, err = mtx.WriteString(strconv.Itoa(geneIndex+1) + "\t" + strconv.Itoa(cellIndex+1) + "\t" + strconv.Itoa(count) + "\n")
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

	for i, geneName := range genes {
		_, err = genesTSV.WriteString(geneIDs[i] + "\t" + geneName + "\n")
		if err != nil {
			log.Fatal(err)
		}
	}
	genesTSV.Flush()
}

// CountMtx converts matrix.tcc.mtx file to gene counts
func CountMtx(tccFile string) {

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
	var geneIDs []string
	var geneIndices []int
	var genes []string
	for _, v := range geneMap {
		geneIDs = append(geneIDs, v[0])
		genes = append(genes, v[1])
	}
	geneIDs, geneIndices = uniqueWithIndices(geneIDs)
	genes = subset(genes, geneIndices)

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
	tccMtx, err := os.Open(tccFile)
	if err != nil {
		log.Fatal(err)
	}
	defer tccMtx.Close()

	// Scan .mtx file
	tccScanner := bufio.NewScanner(tccMtx)

	var line string
	headerRead := false

	var lineIndex int
	var nCounts uint64

	for tccScanner.Scan() {
		line = tccScanner.Text()

		// Skip header and comments
		if line[0:1] == "%" {
			continue
		}
		// Skip the line of the Matrix Market file defining dimensions
		if !headerRead {
			headerRead = true
			continue
		}

		s := strings.Split(line, "\t")

		cell, err := strconv.Atoi(s[0])
		if err != nil {
			log.Fatal(err)
		}
		cell = cell - 1 // MM is 1-based format

		ec, err := strconv.Atoi(s[1])
		if err != nil {
			log.Fatal(err)
		}
		ec = ec - 1 // MM is 1-based format

		if gene, ok := ecGeneMap[int32(ec)]; ok {
			count, err := strconv.Atoi(s[2])
			if err != nil {
				log.Fatal(err)
			}
			// Note: this only works when a single gene corresponds to a single EC
			geneInt := index(genes, gene[0])

			// Consider an option to add up counts
			if calcReadCount {
				cellGeneCounts[cell][geneInt] += uint32(count)
				nCounts += uint64(count)
			} else {
				cellGeneCounts[cell][geneInt]++
				nCounts++
			}
		}
		lineIndex++
		fmt.Printf("\rProcessed line %d, recorded %d counts", lineIndex, nCounts)
	}
	fmt.Printf("\n")

	// Create output dir if does not exist
	if _, err := os.Stat(outputDir); os.IsNotExist(err) {
		os.Mkdir(outputDir, 0700)
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

	_, err = mtx.WriteString(strconv.Itoa(len(genes)) + "\t" + strconv.Itoa(len(cells)) + "\t" + strconv.Itoa(nElements) + "\n")
	if err != nil {
		log.Fatal(err)
	}

	var count int
	var cellCounts []uint32
	for cellIndex := range cells {
		cellCounts = cellGeneCounts[cellIndex]
		for geneIndex := range genes {
			count = int(cellCounts[geneIndex])
			if count > 0 {
				_, err = mtx.WriteString(strconv.Itoa(geneIndex+1) + "\t" + strconv.Itoa(cellIndex+1) + "\t" + strconv.Itoa(count) + "\n")
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

	for i, geneName := range genes {
		_, err = genesTSV.WriteString(geneIDs[i] + "\t" + geneName + "\n")
		if err != nil {
			log.Fatal(err)
		}
	}
	genesTSV.Flush()
}
