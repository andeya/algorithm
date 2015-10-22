// Copyright 2011, 2012  Johann Höchtl. All rights reserved.
// Use of this source code is governed by a Modified BSD License
// that can be found in the LICENSE file.

// This package provides functionality to calculate intermediate
// steps for division "the pen and paper method". Thus it may support
// controlling intermediate steps when pupils start learning to divide.
package schoolcalc

import (
	"fmt"
	"math/big"
	"strings"
)

// Calculation specifiers affecting precision once the dividend is exhausted, but the remainder is not (yet) zero
const (
	SDivPrecReached uint8 = 1 << 7               // a flag; if set, division stops once the remainder is zero
	SDivPrecDefault       = SDivPrecReached | 10 // per default, calculate until the remainder reaches 0 or a maximum of 10 iterations past dividend
	SDivPrecMax           = 127                  // 127 digits maximum precision
)

type Iremainder struct {
	Indent     int
	Iremainder string
}

type SDivide struct {
	Dividend, Divisor,
	NormalizedDividend, NormalizedDivisor,
	Result, Remainder string
	DivisionSteps    []Iremainder
	Prec, ActualPrec uint8
	Exact            bool
	Negative         bool
}

// SchoolDivide accepts two strings and will return the result and intermediate remainders.
//
// dividend and divisor are strings, both may contain fractions denoted by '.'
// prec set's the precision, see Calculation specifiers
func SchoolDivide(dividend, divisor string, prec uint8) (sd *SDivide, err error) {

	var dividendsuffixlen, divisorsuffixlen int
	var endresult string
	var runningprec uint8
	var dividendep int

	mydividend := dividend
	mydivisor := divisor
	first := true
	exact := false
	steps := []Iremainder{}
	onestep := Iremainder{}
	negative := false

	if len(mydividend) <= 0 {
		return nil, fmt.Errorf("Dividend must not be null")
	}

	if len(mydivisor) <= 0 {
		return nil, fmt.Errorf("Divisor must not be null")
	}

	if mydividend[0] == '-' {
		negative = !negative
		mydividend = mydividend[1:]
	}

	if mydivisor[0] == '-' {
		negative = !negative
		mydivisor = mydivisor[1:]
	}

	// check if the divisor contains a fraction
	splitstrings := strings.Split(mydivisor, ".")
	if slen := len(splitstrings); slen > 2 {
		return nil, fmt.Errorf("Not a valid divisor: \"%s\"", divisor)
	} else if slen == 2 {
		// if it contains a suffix, we determine the length of it
		divisorsuffixlen = len(splitstrings[1])
		mydivisor = splitstrings[1]
		// treatment of 0.xyz numerals; if the significant part doesn't start with "0" it can be a part of the numeral
		if splitstrings[0] != "0" {
			mydivisor = splitstrings[0] + mydivisor
		}
	}

	// check if the dividend contains a fraction
	splitstrings = strings.Split(mydividend, ".")
	if slen := len(splitstrings); slen > 2 {
		return nil, fmt.Errorf("Not a valid dividend: \"%s\"", dividend)
	} else if slen == 2 {
		dividendsuffixlen = len(splitstrings[1])
		mydividend = splitstrings[1]
		if splitstrings[0] != "0" {
			mydividend = splitstrings[0] + mydividend
		}
	}

	// How many zeros have to be appended to get rid of fractional numeric parts both for dividend and divisor
	padlen := dividendsuffixlen - divisorsuffixlen

	// depending on padlenght, we have to fill the dividend or the divisor with padlen zeros
	if padlen < 0 {
		mydividend += strings.Repeat("0", -padlen)
	} else if padlen > 0 {
		mydivisor += strings.Repeat("0", padlen)
	}

	bigdivisor, ok := big.NewInt(0).SetString(mydivisor, 10)
	if !ok {
		return nil, fmt.Errorf("Not a divisor: \"%s\"", mydivisor)
	}

	bigintermediatedividend := big.NewInt(0)
	onebig := big.NewInt(0)

	// start to divide as soon as the dividend is greater than the divisor
	for dividendep = 1; ; dividendep++ {
		if _, ok = bigintermediatedividend.SetString(mydividend[0:dividendep], 10); !ok {
			return nil, fmt.Errorf("Not a dividend: \"%s\"", dividend)
		}

		// We are done checking once the running dividend is bigger than the divisor or if it can't get bigger because dividend < divisor
		if !(bigintermediatedividend.Cmp(bigdivisor) < 0 && dividendep < len(mydividend)) {
			break
		}
	}

	for {
		intermediateresult := big.NewInt(0).Div(bigintermediatedividend, bigdivisor)
		endresult += intermediateresult.String()
		bigintermediatedividend = big.NewInt(0).Mul(big.NewInt(10), big.NewInt(0).Rem(bigintermediatedividend, bigdivisor))

		if dividendep < len(mydividend) {
			// if the dividend is not exhausted, we have to add the next position of the dividend to the running dividend
			if _, ok = onebig.SetString(string(mydividend[dividendep]), 10); !ok {
				return nil, fmt.Errorf("Not a number: \"%c\" (in %s)", mydividend[dividendep], mydividend)
			}

			if bigintermediatedividend.Cmp(big.NewInt(0)) == 0 {
				onestep.Iremainder = "0"
			} else {
				onestep.Iremainder = ""
			}

			bigintermediatedividend.Add(bigintermediatedividend, onebig)

			// the running dividend is an intermediate step to record
			onestep.Iremainder += bigintermediatedividend.String()
			onestep.Indent = dividendep - len(onestep.Iremainder) + 1
			steps = append(steps, onestep)

			dividendep++

		} else if (SDivPrecMax&prec)-runningprec > 0 {
			// if the dividend is exhausted, calculation continues ...
			onestep.Iremainder = bigintermediatedividend.String()

			// ... until we reach the maximum desired precision or the remainder(= running dividend) is zero
			if bigintermediatedividend.Cmp(big.NewInt(0)) == 0 && (prec&SDivPrecReached) > 0 {
				exact = true

				onestep.Indent = dividendep + int(runningprec) - len(onestep.Iremainder)
				steps = append(steps, onestep)

				break
			} else {
				if onestep.Iremainder == "0" {
					onestep.Iremainder = "00"
				}
			}
			onestep.Indent = dividendep + int(runningprec) - len(onestep.Iremainder) + 1
			steps = append(steps, onestep)

			// the first time we exhaust the running divided, the result will be fractional
			if first {
				endresult += "."
				first = false
			}
			runningprec++
		} else {
			onestep.Iremainder = bigintermediatedividend.Div(bigintermediatedividend, big.NewInt(10)).String()
			onestep.Indent = dividendep + int(runningprec) - len(onestep.Iremainder)
			steps = append(steps, onestep)
			break
		}
	}

	if negative {
		endresult = "-" + endresult
	}

	return &SDivide{Dividend: dividend, Divisor: divisor, Result: endresult, Remainder: onestep.Iremainder,
		NormalizedDividend: mydividend, NormalizedDivisor: mydivisor,
		DivisionSteps: steps,
		Prec:          prec, ActualPrec: runningprec,
		Exact:    exact,
		Negative: negative}, nil
}

func inputSDivisequaltoResultSDiv(inputdivisor, inputdividend, outputdivisor, outputdivided string) bool {
	return inputdividend == outputdivided && inputdivisor == outputdivisor
}

func (sd *SDivide) String() string {
	var blank string
	sresult := fmt.Sprintf("%s : %s = %s\n", sd.Dividend, sd.Divisor, sd.Result)
	if !inputSDivisequaltoResultSDiv(sd.Dividend, sd.Divisor, sd.NormalizedDividend, sd.NormalizedDivisor) {
		sresult += fmt.Sprintf("%s : %s = %s\n", sd.NormalizedDividend, sd.NormalizedDivisor, sd.Result)
	}
	for _, elm := range sd.DivisionSteps {
		blank = strings.Repeat(" ", elm.Indent)
		sresult += fmt.Sprintf("%s%s\n", blank, elm.Iremainder)
	}
	return sresult
}

// A 'Zapfen' consists of 8 multiplications and 8 divisions.
// Example:
//
//         27 * 2 = 54
//         54 * 3 = 162
//        162 * 4 = 648
//        648 * 5 = 3240
//       3240 * 6 = 19440
//      19440 * 7 = 136080
//     136080 * 8 = 1088640
//    1088640 * 9 = 9797760
//    9797760 / 2 = 4898880
//    4898880 / 3 = 1632960
//    1632960 / 4 = 408240
//     408240 / 5 = 81648
//      81648 / 6 = 13608
//      13608 / 7 = 1944
//       1944 / 8 = 243
//        243 / 9 = 27
//

// The struct stores the eight intermediary multiplications, the eight divisions and
// the string length of the longest product to allow proper result formatting.
type Zapfen struct {
	Zapfenzahl *big.Int
	Multzapfen [8]*big.Int
	Divzapfen  [8]*big.Int
	Longest    int
}

func ZapfenRechnung(zapfenzahl *big.Int) (rv *Zapfen) {
	rv = &Zapfen{Zapfenzahl: zapfenzahl}

	// eight multiplications, starting with "2"
	for i := 2; i < 10; i++ {
		rv.Multzapfen[i-2] = new(big.Int)
		rv.Multzapfen[i-2].Mul(zapfenzahl, big.NewInt(int64(i)))
		zapfenzahl = rv.Multzapfen[i-2]
	}

	// calculate the string length of the longest product to allow proper result formatting
	rv.Longest = len(rv.Multzapfen[7].String())
	if rv.Longest == 0 {
		rv.Longest = 1
	}

	// perform eight divisions, starting to divide the last product of the preceding calculation with "2".
	// The end result will be calling zapfenzahl
	for i := 2; i < 10; i++ {
		rv.Divzapfen[i-2] = new(big.Int)
		rv.Divzapfen[i-2].Div(zapfenzahl, big.NewInt(int64(i)))
		zapfenzahl = rv.Divzapfen[i-2]
	}
	return
}

func (rv *Zapfen) String() string {
	var sresult string
	input := rv.Zapfenzahl
	for i := 0; i < 8; i++ {
		sresult += fmt.Sprintf("%*d * %d = %d\n", rv.Longest, input, i+2, rv.Multzapfen[i])
		input = rv.Multzapfen[i]
	}

	for i := 0; i < 8; i++ {
		sresult += fmt.Sprintf("%*d / %d = %d\n", rv.Longest, input, i+2, rv.Divzapfen[i])
		input = rv.Divzapfen[i]

	}
	return sresult
}
