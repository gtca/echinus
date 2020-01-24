package cmd

import (
	"fmt"

	"github.com/spf13/cobra"
)

func init() {
	rootCmd.AddCommand(versionCmd)
}

var versionCmd = &cobra.Command{
	Use:   "version",
	Short: "Print the version number of echinus",
	Long:  "This command is to be used to get to know the version of echinus",
	Run: func(cmd *cobra.Command, args []string) {
		fmt.Println("echinus v0.1.0")
	},
}
