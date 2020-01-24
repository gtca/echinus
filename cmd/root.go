package cmd

import (
	"fmt"
	"os"

	"github.com/spf13/cobra"
)

var rootCmd = &cobra.Command{
	Use:   "echinus",
	Short: "echinus converts transcript compatibility counts to gene counts",
	Long:  "Echinus provides a fast and easy way to summarise equivalence counts to transcript or gene level.\nWritten in Go.",
}

// Execute the root command
func Execute() {
	if err := rootCmd.Execute(); err != nil {
		fmt.Println(err)
		os.Exit(1)
	}
}
