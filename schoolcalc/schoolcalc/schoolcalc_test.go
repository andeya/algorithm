// Copyright 2011 Johann Höchtl. All rights reserved.
// Use of this source code is governed by a Modified BSD License
// that can be found in the LICENSE file.

// Automated tests for the sdivision package
package schoolcalc

import (
	"testing"
)

// ## TestSchoolDivision
type schoolDivisionParam struct {
	dividend, divisor string
	prec              uint8
}

type schoolDivisionTest struct {
	in  schoolDivisionParam
	out *SDivide
}

var schoolDivisionTests = []schoolDivisionTest{
	{schoolDivisionParam{"100", "5", SDivPrecReached | 8},
		&SDivide{Dividend: "100", Divisor: "5", NormalizedDividend: "100", NormalizedDivisor: "5", Result: "20", Remainder: "0", DivisionSteps: []Iremainder{Iremainder{Indent: 1, Iremainder: "00"}, Iremainder{Indent: 2, Iremainder: "0"}}, Prec: 0x88, Exact: true, Negative: false}},
	{schoolDivisionParam{"100.5", "5", SDivPrecReached | 8},
		&SDivide{Dividend: "100.5", Divisor: "5", NormalizedDividend: "1005", NormalizedDivisor: "50", Result: "20.1", Remainder: "0", DivisionSteps: []Iremainder{Iremainder{Indent: 2, Iremainder: "05"}, Iremainder{Indent: 3, Iremainder: "50"}, Iremainder{Indent: 4, Iremainder: "0"}}, Prec: 0x88, Exact: true, Negative: false}},
	{schoolDivisionParam{"100.5", "5.5", SDivPrecReached | 8},
		&SDivide{Dividend: "100.5", Divisor: "5.5", NormalizedDividend: "1005", NormalizedDivisor: "55", Result: "18.27272727", Remainder: "15", DivisionSteps: []Iremainder{Iremainder{Indent: 1, Iremainder: "455"}, Iremainder{Indent: 2, Iremainder: "150"}, Iremainder{Indent: 3, Iremainder: "400"}, Iremainder{Indent: 4, Iremainder: "150"}, Iremainder{Indent: 5, Iremainder: "400"}, Iremainder{Indent: 6, Iremainder: "150"}, Iremainder{Indent: 7, Iremainder: "400"}, Iremainder{Indent: 8, Iremainder: "150"}, Iremainder{Indent: 9, Iremainder: "400"},
			Iremainder{Indent: 10, Iremainder: "15"}}, Prec: 0x88, Exact: false, Negative: false}},
	{schoolDivisionParam{"-100.5", "5.56", SDivPrecReached | 2},
		&SDivide{Dividend: "-100.5", Divisor: "5.56", NormalizedDividend: "10050", NormalizedDivisor: "556", Result: "-18.07", Remainder: "308", DivisionSteps: []Iremainder{Iremainder{Indent: 1, Iremainder: "4490"}, Iremainder{Indent: 3, Iremainder: "420"}, Iremainder{Indent: 3, Iremainder: "4200"}, Iremainder{Indent: 4, Iremainder: "308"}}, Prec: 0x82, Exact: false, Negative: true}},
	{schoolDivisionParam{"5", "100", SDivPrecReached | 2},
		&SDivide{Dividend: "5", Divisor: "100", NormalizedDividend: "5", NormalizedDivisor: "100", Result: "0.05", Remainder: "0", DivisionSteps: []Iremainder{Iremainder{Indent: 0, Iremainder: "50"}, Iremainder{Indent: 0, Iremainder: "500"}, Iremainder{Indent: 2, Iremainder: "0"}}, Prec: 0x82, Exact: false, Negative: false}},
	{schoolDivisionParam{"-5", "100", SDivPrecReached | 2},
		&SDivide{Dividend: "-5", Divisor: "100", NormalizedDividend: "5", NormalizedDivisor: "100", Result: "-0.05", Remainder: "0", DivisionSteps: []Iremainder{Iremainder{Indent: 0, Iremainder: "50"}, Iremainder{Indent: 0, Iremainder: "500"}, Iremainder{Indent: 2, Iremainder: "0"}}, Prec: 0x82, Exact: false, Negative: true}},
	{schoolDivisionParam{"2", "0.5", SDivPrecReached | 2},
		&SDivide{Dividend: "2", Divisor: "0.5", NormalizedDividend: "20", NormalizedDivisor: "5", Result: "4", Remainder: "0", DivisionSteps: []Iremainder{Iremainder{Indent: 1, Iremainder: "0"}}, Prec: 0x82, Exact: true, Negative: false}},
	{schoolDivisionParam{"0.5", "0.5", SDivPrecReached | 2},
		&SDivide{Dividend: "0.5", Divisor: "0.5", NormalizedDividend: "5", NormalizedDivisor: "5", Result: "1", Remainder: "0", DivisionSteps: []Iremainder{Iremainder{Indent: 0, Iremainder: "0"}}, Prec: 0x82, Exact: true, Negative: false}},
	{schoolDivisionParam{"10065767", "55.7", SDivPrecReached | 2},
		&SDivide{Dividend: "10065767", Divisor: "55.7", NormalizedDividend: "100657670", NormalizedDivisor: "557", Result: "180713.94", Remainder: "542", DivisionSteps: []Iremainder{Iremainder{Indent: 1, Iremainder: "4495"}, Iremainder{Indent: 3, Iremainder: "397"}, Iremainder{Indent: 3, Iremainder: "3976"}, Iremainder{Indent: 5, Iremainder: "777"}, Iremainder{Indent: 5, Iremainder: "2200"}, Iremainder{Indent: 6, Iremainder: "5290"}, Iremainder{Indent: 7, Iremainder: "2770"}, Iremainder{Indent: 8, Iremainder: "542"}}, Prec: 0x82, Exact: false, Negative: false}},
	{schoolDivisionParam{"100X65767", "55.7", SDivPrecReached | 2},
		nil},
	{schoolDivisionParam{"100X65767", "Y.7", SDivPrecReached | 2},
		nil},
	{schoolDivisionParam{"", "Y.7", SDivPrecReached | 2},
		nil},
	{schoolDivisionParam{"12.6546.65", "Y.7", SDivPrecReached | 2},
		nil},
}

func devidecompare(p1, p2 *SDivide) bool {

	if p1.NormalizedDividend != p2.NormalizedDividend {
		return false
	}
	if p1.NormalizedDivisor != p2.NormalizedDivisor {
		return false
	}
	if p1.Result != p2.Result {
		return false
	}
	if p1.Remainder != p2.Remainder {
		return false
	}
	if p1.Exact != p2.Exact {
		return false

	}
	if p1.Negative != p2.Negative {
		return false
	}

	if len(p1.DivisionSteps) != len(p2.DivisionSteps) {
		return false
	}

	for index, elm1 := range p1.DivisionSteps {
		if elm1.Indent != p2.DivisionSteps[index].Indent {
			return false
		}
		if elm1.Iremainder != p2.DivisionSteps[index].Iremainder {
			return false
		}
	}
	return true
}

func TestSchoolDivision(t *testing.T) {
	for cnt, test := range schoolDivisionTests {
		out, err := SchoolDivide(test.in.dividend, test.in.divisor, test.in.prec)
		if err == nil && test.out != nil {
			if !devidecompare(out, test.out) {
				t.Errorf("SchoolDivide [%d]: Expected %v, got %v", cnt, test.out, out)
			}
		} else if err != nil && test.out != nil {
			t.Errorf("SchoolDivide [%d]: Got error %s but test case says it's save to divide %v", cnt, err.Error(), test.in)
		} else if err == nil && test.out == nil {
			t.Errorf("SchoolDivide [%d]: Got a result %v while test case says it's erroneous to divide %v", cnt, test.out, test.in)
		}
	}
}
