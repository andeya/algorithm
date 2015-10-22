// Copyright (c) Harri Rautila, 2012

// This file is part of go.opt/matrix package. It is free software, distributed
// under the terms of GNU Lesser General Public License Version 3, or any later
// version. See the COPYING tile included in this archive.

package matrix

import (
    //"errors"
    "fmt"
    "strconv"
    "strings"
)

// Convert matrix to row-major string representation. 
func (A *FloatMatrix) String() string {
    s := ""
    if A == nil {
        return "<nil>"
    }
    step := A.LeadingIndex()
    for i := 0; i < A.Rows(); i++ {
        if i > 0 {
            s += "\n"
        }
        s += "["
        for j := 0; j < A.Cols(); j++ {
            if j > 0 {
                s += ", "
            }
            s += fmt.Sprintf("%9.2e", A.elements[j*step+i])
        }
        s += "]"
    }
    return s
}

// Convert matrix to row-major string representation using format as element format. 
func (A *FloatMatrix) ToString(format string) string {
    s := ""
    if A == nil {
        return "<nil>"
    }
    step := A.LeadingIndex()
    for i := 0; i < A.Rows(); i++ {
        if i > 0 {
            s += "\n"
        }
        s += "["
        for j := 0; j < A.Cols(); j++ {
            if j > 0 {
                s += ", "
            }
            s += fmt.Sprintf(format, A.elements[j*step+i])
        }
        s += "]"
    }
    return s
}

// Convert matrix to row-major string representation. 
func (A *FloatMatrix) ConvertToString() string {
    s := ""
    if A == nil {
        return "<nil>"
    }
    step := A.LeadingIndex()
    for i := 0; i < A.Rows(); i++ {
        if i > 0 {
            s += "\n"
        }
        s += "["
        for j := 0; j < A.Cols(); j++ {
            if j > 0 {
                s += ", "
            }
            s += fmt.Sprintf("%.17f", A.elements[j*step+i])
        }
        s += "]"
    }
    return s
}

func FloatParse(s string) (A *FloatMatrix, err error) {
    if strings.Index(s, ";") != -1 {
        // Matlab style [a b c; d e f] recognised from semicolon
        return floatParseMatLab(s)
    } else if strings.Index(s, "{") != -1 {
        // special format {rows cols [data]}
        return floatParseSpe(s)
    }
    // finaly try Python string [row]\n[row]
    return floatParsePy(s)
}

// Parse a matlab-style row-major matrix representation eg [a b c; d e f]
// and return a DenseFLoatMatrix.
func floatParseMatLab(s string) (A *FloatMatrix, err error) {
    var arrays [][]float64
    start := strings.Index(s, "[")
    end := strings.LastIndex(s, "]")
    if start == -1 || end == -1 {
        err = onError("Unrecognized matrix string")
        return
    }
    rowStrings := strings.Split(s[start+1:end], ";")
    //nrows := len(rowStrings)
    ncols := 0
    for _, row := range rowStrings {
        rowElems := strings.Split(strings.Trim(row, " "), " ")
        if ncols == 0 {
            ncols = len(rowElems)
        } else if ncols != len(rowElems) {
            err = ErrorDimensionMismatch
            return
        }
        row := []float64{}
        for _, valString := range rowElems {
            var val float64
            val, err = strconv.ParseFloat(valString, 64)
            if err != nil {
                return
            }
            row = append(row, val)
        }
        arrays = append(arrays, row)
    }
    A = FloatMatrixFromTable(arrays, RowOrder)
    return
}

// Parse python cvxopt string representation of a matrix.
//   [1,0 2.0 3.0]
//   [1.1 2.1 3.1]
// Returns a new FloatMatrix.
func floatParsePy(s string) (A *FloatMatrix, err error) {
    var arrays [][]float64
    // rowString is matrix row starting with '[' character.
    // Remove newlines and split on ']'
    rowStrings := strings.Split(strings.Trim(s, "\n"), "]")
    //fmt.Printf("rows string: '%v'\n", rowStrings)
    ncols := 0
    firstRow := true
    currow := 0
    hasComma := strings.Index(s, ",") != -1
rows:
    for _, rowStr := range rowStrings {
        if len(rowStr) == 0 {
            continue rows
        }
        rowStr := strings.Trim(rowStr, " \n][")
        if hasComma {
            // replace commas with ws
            rowStr = strings.Replace(rowStr, ",", " ", -1)
        }
        rowElems := strings.Fields(rowStr)
        row := []float64{}
        collen := 0
        for _, valString := range rowElems {
            var val float64
            if firstRow {
                ncols += 1
            }
            collen += 1
            val, err = strconv.ParseFloat(valString, 64)
            if err != nil {
                return
            }
            row = append(row, val)
        }
        currow += 1
        if firstRow {
            firstRow = false
        }
        if collen != ncols {
            err = onError(fmt.Sprintf("row %d: num columns %d, expected %d\n", currow, collen, ncols))
            return
        }
        arrays = append(arrays, row)
    }
    A = FloatMatrixFromTable(arrays, RowOrder)
    return
}

// Parse string in format '{nrows ncols [element list]}'
func floatParseSpe(s string) (A *FloatMatrix, err error) {
    err = nil
    A = nil

    start := strings.Index(s, "{")
    end := strings.LastIndex(s, "}")
    if start == -1 || end == -1 {
        err = onError("Unrecognized matrix string")
        return
    }
    dstart := strings.Index(s, "[")
    dend := strings.Index(s, "]")
    if dstart == -1 || dend == -1 {
        err = onError("Unrecognized matrix string")
        return
    }

    var nrows, ncols int64
    var val float64
    sizes := strings.Fields(strings.Replace(s[start+1:dstart], ",", " ", -1))
    nrows, err = strconv.ParseInt(sizes[0], 10, 64)
    if err != nil {
        return
    }
    ncols, err = strconv.ParseInt(sizes[1], 10, 64)
    if err != nil {
        return
    }
    if dend-dstart == 1 {
        // empty matrix; no rows
        A = FloatZeros(0, 1)
        return
    }
    es := s[dstart+1 : dend]
    if strings.Index(es, ",") != -1 {
        es = strings.Replace(es, ",", " ", -1)
    }
    elems := strings.Fields(es)
    data := make([]float64, len(elems))
    for k, elem := range elems {
        val, err = strconv.ParseFloat(strings.Trim(elem, " "), 64)
        if err != nil {
            return
        }
        data[k] = val
    }
    A = FloatNew(int(nrows), int(ncols), data)
    return
}

// Local Variables:
// tab-width: 4
// End:
