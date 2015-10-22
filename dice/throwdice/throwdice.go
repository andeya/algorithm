package main

import (
	"bufio"
	"flag"
	"fmt"
	"github.com/henrylee2cn/algorithm/dice"
	"io"
	"os"
	"strconv"
	"strings"
)

func intsToStrings(ints []int) (strings []string) {
	strings = make([]string, len(ints))
	for i, v := range ints {
		strings[i] = strconv.Itoa(v)
	}
	return
}

func usage() {
	fmt.Fprintf(os.Stderr, "Usage of %s:\n", os.Args[0])
	fmt.Fprintf(os.Stderr, "  "+os.Args[0]+" [<roll description>...]\n")
	fmt.Fprintf(os.Stderr, "  -help\tprints this help message\n")
	flag.PrintDefaults()
}

var d6 bool

func init() {
	flag.Usage = usage
	flag.BoolVar(&d6, "d6", false, "use the d6 system")
}

func main() {
	flag.Parse()

	if flag.NArg() < 1 {
		reader := bufio.NewReader(os.Stdin)
		for {
			line := ""
			if buf, pre, err := reader.ReadLine(); err != nil {
				if err != io.EOF {
					fmt.Println(err)
				} else {
					break
				}
			} else {
				line = line + string(buf)
				if !pre {
					lowerLine := strings.ToLower(line)
					if lowerLine == "exit" || lowerLine == "quit" {
						return
					}
					printDiceRoll(line)
					line = ""
				}
			}
		}
	} else {
		for i := 0; i < flag.NArg(); i++ {
			printDiceRoll(flag.Arg(i))
		}
	}
}

func printDiceRoll(description string) {
	if d6 {
		if number, err := strconv.Atoi(description); err != nil {
			fmt.Println(err)
		} else {
			roll := dice.RollD6(number)
			fmt.Printf("%v (%v)\n", roll.Total, strings.Join(intsToStrings(roll.Rolls), ", "))
		}
	} else {
		if roll, err := dice.Roll(description); err != nil {
			fmt.Println(err)
		} else {
			fmt.Printf("%v (%v)\n", roll.Total, strings.Join(intsToStrings(roll.Rolls), ", "))
		}
	}
}
