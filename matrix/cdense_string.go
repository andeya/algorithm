// Copyright (c) Harri Rautila, 2012

// This file is part of go.opt/matrix package. It is free software, distributed
// under the terms of GNU Lesser General Public License Version 3, or any later
// version. See the COPYING tile included in this archive.

package matrix

import (
    "errors"
    "fmt"
    "strconv"
    "strings"
)

// Convert matrix to row-major string representation. 
func (A *ComplexMatrix) String() string {
    s := ""
    step := A.LeadingIndex()
    for i := 0; i < A.Rows(); i++ {
        if i > 0 {
            s += "\n"
        }
        s += "["
        for j := 0; j < A.Cols(); j++ {
            if j > 0 {
                s += " "
            }
            s += fmt.Sprintf("%v", A.elements[j*step+i])
        }
        s += "]"
    }
    return s
}

// Parse a matlab-style row-major matrix representation eg [a b c; d e f]
// and return a Matrix. Each element is pair of floats e.g.
// [(1.0+0i) (0.0+0i); (4-2i) (-1-2i)]
func ComplexParse(s string) (A *ComplexMatrix, err error) {
    var arrays [][]complex128
    start := strings.Index(s, "[")
    end := strings.LastIndex(s, "]")
    if start == -1 || end == -1 {
        err = errors.New("Unrecognized matrix string")
        return
    }
    rowStrings := strings.Split(s[start+1:end], ";")
    //nrows := len(rowStrings)
    ncols := 0
    for _, row := range rowStrings {
        rowElems := strings.Split(strings.Trim(row, ") "), ")")
        if ncols == 0 {
            ncols = len(rowElems)
        } else if ncols != len(rowElems) {
            err = ErrorDimensionMismatch
            return
        }
        row := []complex128{}
        for _, valString := range rowElems {
            var cval complex128
            valString = strings.Trim(valString, "( ")
            cval, err = parseComplex(valString)
            if err != nil {
                return
            }
            row = append(row, cval)
        }
        arrays = append(arrays, row)
    }
    A = ComplexMatrixFromTable(arrays, RowOrder)
    return
}

func parseComplex(s string) (v complex128, err error) {
    var re, im string
    var rval, ival float64

    v = complex(0, 0)
    str := strings.Trim(s, " ")
    ri := 0
    if str[0] == '-' || str[0] == '+' {
        ri += 1
    }
    if n := strings.Index(str[ri:], "+"); n != -1 {
        re = strings.Trim(str[:n+ri], " ")
        im = strings.Trim(str[n+ri:], " ")
    } else if n := strings.Index(str[ri:], "-"); n != -1 {
        re = strings.Trim(str[:n+ri], " ")
        im = strings.Trim(str[n+ri:], " ")
    } else {
        err = errors.New("invalid complex number")
        return
    }
    //fmt.Printf("re: '%s', im: '%s'\n", re, im)
    if im[len(im)-1] != 'i' {
        err = errors.New("invalid complex number, no imaginary part")
        return
    }
    rval, err = strconv.ParseFloat(re, 64)
    if err != nil {
        return
    }
    ival, err = strconv.ParseFloat(im[:len(im)-1], 64)
    if err != nil {
        return
    }
    v = complex(rval, ival)
    return
}

// Local Variables:
// tab-width: 4
// End:
