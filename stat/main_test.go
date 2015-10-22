package stat

/* main_test.go
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Jim Davies, Brian Gough
 * Copyright (C) 2012, 2013 G.vd.Schoot
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

import (
	"fmt"
	"math"
	"os"
	"testing"
)

func gsl_test(status bool, test_description string, a ...interface{}) {

	if status == true {
		s := fmt.Sprintf(test_description, a...)
		fmt.Printf("FAIL: %s\n", s)
		os.Exit(1)
	}
}

func gsl_test_rel(result, expected, relative_error float64, test_description string, a ...interface{}) {

	var status bool

	// Check for NaN vs inf vs number

	if math.IsNaN(result) || math.IsNaN(expected) {
		status = math.IsNaN(result) != math.IsNaN(expected)
	} else if math.IsInf(result, 0) || math.IsInf(expected, 0) {
		status = math.IsInf(result, 0) != math.IsInf(expected, 0)
	} else if expected != 0 {
		status = math.Abs(result-expected)/math.Abs(expected) > relative_error
	} else {
		status = math.Abs(result) > relative_error
	}

	if status == true {
		s := fmt.Sprintf(test_description, a...)
		fmt.Printf("FAIL: %s (%f observed vs %g expected)\n",
			s, result, expected)
		os.Exit(1)
	}
}

// Main test routine
func TestMain(t *testing.T) {

	for s1 := 1; s1 < 4; s1++ {
		s2 := 1
		if s1 >= 3 {
			s2 = s1 - 1
		}

		testFloat64Slice(t, s1, s2)
		testIntSlice(t, s1, s2)
	}

	testNist()
}
