package cmd

import (
	"testing"
)

func TestBinaryToString(t *testing.T) {
	want := "AAAAAA"
	if got := binaryToString(0, uint32(len(want))); got != want {
		t.Errorf("0 => %q, want %q", got, want)
	}

	want = "AGCTAG"
	if got := binaryToString(626, uint32(len(want))); got != want {
		t.Errorf("626 => %q, want %q", got, want)
	}
}
