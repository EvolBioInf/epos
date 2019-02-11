package main

import (
	"os"
	"strings"

	"github.com/evolbioinf/epos2plot/lib"
)

func run(f *os.File, opts lib.Opts) {
	if opts.R {
		lib.PrintRaw(f)
	} else {
		data := lib.Read(f)
		quant := lib.Quantiles(data, opts)
		lib.PrintQuantiles(quant)
	}
}
var VERSION, DATE string
func main() {
	date := strings.Replace(DATE, "_", " ", -1)
	opts := lib.ParseCL(date, VERSION)
	if len(opts.Files) == 0 {
		run(os.Stdin, opts)
	} else {
		for _, f := range opts.Files {
			fi, err := os.Open(f)
			lib.Check(err)
			run(fi, opts)
			err = fi.Close()
			lib.Check(err)
		}
	}
}
