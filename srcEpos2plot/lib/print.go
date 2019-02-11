package lib

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"strings"
	"strconv"
)

func PrintQuantiles(qu []Quantile) {
	fmt.Printf("#Time\tLowerQ\tMedian\tUpperQ\n")
	for i, q := range qu {
		if q.T < math.Inf(1) && (i == 0 || qu[i].T != qu[i-1].T) {
			fmt.Printf("%g\t%g\t%g\t%g\n", q.T, q.L, q.M, q.U)
		} 
	}
}

func PrintRaw(f *os.File) {
	var st, et float64  // Start and end times for interval.
	var pn, cn float64 // previous and current population size
	var er error
	input := bufio.NewScanner(f)
	first := true
	veryFirst := true
	for input.Scan() {
		s := input.Text()
		if s[0:1] == "#" {    // Reset start time.
			st = 0
			if !first {
				fmt.Printf("%g\t%g\n", et, cn)
			}
			first = true
		} else {
			arr := strings.Split(input.Text(), "\t")
			et, er = strconv.ParseFloat(arr[1], 64) // Get time for end of interval.
			Check(er)
			cn, er = strconv.ParseFloat(arr[2], 64) // Get population size.
			Check(er)
			if st > 0 && et > 0 {
				if first {
					first = false
					if veryFirst {
						fmt.Printf("%g\t%g\n", 0., pn)
						veryFirst = false
					} else {
						fmt.Printf("\n%g\t%g\n", 0., pn)
					}
				}
				fmt.Printf("%g\t%g\n", st, pn)
				fmt.Printf("%g\t%g\n", st, cn)
			}
			st = et
			pn = cn
		}
	}
	fmt.Printf("%g\t%g\n", et, cn)
}
