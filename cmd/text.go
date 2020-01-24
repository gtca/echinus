package cmd

import (
	"github.com/spf13/cobra"
)

// var busFile string

func init() {
	textCmd.Flags().StringVarP(&busFile, "bus", "b", "", "Input BUS file")
	rootCmd.AddCommand(textCmd)
}

var textCmd = &cobra.Command{
	Use:   "text",
	Short: "Get the plain text version of a bustile",
	Long:  "Convert a busfile to a plain text barcode-UMI-EC-count file",
	Run: func(cmd *cobra.Command, args []string) {
		// File with transcript compatibility counts (matrix.tcc.mtx)
		ReadBUS(busFile)
	},
}
