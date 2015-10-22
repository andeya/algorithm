// Copyright (c) Harri Rautila, 2012,2013

// This file is part of go.opt/matrix package. It is free software, distributed
// under the terms of GNU Lesser General Public License Version 3, or any later
// version. See the COPYING tile included in this archive.

package matrix

import (
    "errors"
    "fmt"
)

// Function that are depreceated and will be removed for good.

// Get values for indexed elements. **DEPRECEATED**
func (A *FloatMatrix) GetIndexesFromArray(indexes []int) []float64 {
    return A.GetIndexes(indexes...)
}

// Get copy of i'th row. Row elements are copied to vals array. 
// Returns the array. If vals array is too small new slice is allocated and 
// returned with row elements. **DEPRECEATED**
func (A *FloatMatrix) GetRowArray(i int, vals []float64) []float64 {
    if vals == nil || cap(vals) < A.Cols() {
        vals = make([]float64, A.Cols())
    }
    step := A.LeadingIndex()
    if i < 0 {
        i += A.Rows()
    }
    for j := 0; j < A.Cols(); j++ {
        vals[j] = A.elements[j*step+i]
    }
    return vals
}

// Get copy of i'th row. Return parameter matrix. If vec is too small 
// reallocate new vector and return it. **DEPRECEATED**
// Use SubMatrix function instead.
func (A *FloatMatrix) GetRow(i int, vec *FloatMatrix) *FloatMatrix {
    if vec == nil || vec.NumElements() < A.Cols() {
        vec = FloatZeros(1, A.Cols())
    }
    step := A.LeadingIndex()
    ar := vec.FloatArray()
    if i < 0 {
        i += A.Rows()
    }
    for j := 0; j < A.Cols(); j++ {
        ar[j] = A.elements[j*step+i]
    }
    return vec
}

// Get copy of i'th column. See GetRow. **DEPRECEATED**
// Use SubMatrix function instead.
func (A *FloatMatrix) GetColumn(i int, vec *FloatMatrix) *FloatMatrix {
    if vec == nil || vec.NumElements() < A.Rows() {
        vec = FloatZeros(A.Rows(), 1)
    }
    step := A.LeadingIndex()
    ar := vec.FloatArray()
    if i < 0 {
        i += A.Cols()
    }
    for j := 0; j < A.Rows(); j++ {
        ar[j] = A.elements[i*step+j]
    }
    return vec
}

// Get copy of i'th column. See GetRow. **DEPRECEATED**
func (A *FloatMatrix) GetColumnArray(i int, vec []float64) []float64 {
    if cap(vec) < A.Rows() {
        vec = make([]float64, A.Rows())
    }
    step := A.LeadingIndex()
    if i < 0 {
        i += A.Cols()
    }
    for j := 0; j < A.Rows(); j++ {
        vec[j] = A.elements[i*step+j]
    }
    return vec
}

// Set values of i'th row. **DEPRECEATED**
func (A *FloatMatrix) SetRowArray(i int, vals []float64) {
    step := A.LeadingIndex()
    if i < 0 {
        i = A.Rows() + i
    }
    for j := 0; j < A.Cols(); j++ {
        A.elements[j*step+i] = vals[j]
    }
}

// Set values on i'th row  of columns pointed with cols array. 
// For all j in indexes: A[i,j] = vals[k] where k is j's index in indexes array.
// **DEPRECEATED**
func (A *FloatMatrix) SetAtRowArray(i int, cols []int, vals []float64) {
    step := A.LeadingIndex()
    if i < 0 {
        i = A.Rows() + i
    }
    for k, j := range cols {
        if j < 0 {
            j += A.Cols()
        }
        A.elements[j*step+i] = vals[k]
    }
}

// Set values of i'th row. Matrix vals is either (A.Cols(), 1) or (1, A.Cols()) matrix.
// **DEPRECEATED**
func (A *FloatMatrix) SetRow(i int, vals *FloatMatrix) {
    step := A.LeadingIndex()
    if i < 0 {
        i = A.Rows() + i
    }
    for j := 0; j < A.Cols(); j++ {
        A.elements[j*step+i] = vals.elements[j]
    }
}

// Set values  on i'th row of columns pointed with cols array. 
// For all j in indexes: A[i,j] = vals[j]. Matrix vals is either (A.Cols(),1) or
// (1, A.Cols()) matrix.
// **DEPRECEATED**
func (A *FloatMatrix) SetAtRow(i int, cols []int, vals *FloatMatrix) {
    step := A.LeadingIndex()
    if i < 0 {
        i = A.Rows() + i
    }
    for _, j := range cols {
        if j < 0 {
            j += A.Cols()
        }
        A.elements[j*step+i] = vals.elements[j]
    }
}

// Set values of i'th column. **DEPRECEATED**
func (A *FloatMatrix) SetColumnArray(i int, vals []float64) {
    step := A.LeadingIndex()
    if i < 0 {
        i = A.Cols() + i
    }
    for j := 0; j < A.Rows(); j++ {
        A.elements[i*step+j] = vals[j]
    }
}

// Set values on i'th column of rows pointed by rows array. It assumes
// that len(rows) <= len(vals). **DEPRECEATED**
func (A *FloatMatrix) SetAtColumnArray(i int, rows []int, vals []float64) {
    step := A.LeadingIndex()
    if i < 0 {
        i = A.Cols() + i
    }
    for k, j := range rows {
        if j < 0 {
            j += A.Rows()
        }
        A.elements[i*step+j] = vals[k]
    }
}

// Set values of i'th column. Matrix vals is either (A.Rows(), 1) or (1, A.Rows()) matrix.
// **DEPRECEATED**
func (A *FloatMatrix) SetColumn(i int, vals *FloatMatrix) {
    step := A.LeadingIndex()
    if i < 0 {
        i = A.Cols() + i
    }
    for j := 0; j < A.Rows(); j++ {
        A.elements[i*step+j] = vals.elements[j]
    }
}

// Set values on i'th column of rows pointer by rows array. It assumes
// that  max(rows) < vals.NumElements(). 
// **DEPRECEATED**
func (A *FloatMatrix) SetAtColumn(i int, rows []int, vals *FloatMatrix) {
    step := A.LeadingIndex()
    if i < 0 {
        i = A.Cols() + i
    }
    for _, j := range rows {
        if j < 0 {
            j += A.Rows()
        }
        A.elements[i*step+j] = vals.elements[j]
    }
}

// Set values for sub-matrix starting at (row, col). If row+mat.Rows() greater than
// A.Rows() or col+mat.Cols() greater than A.Cols() matrix A is not changed.
// **DEPRECEATED** (Use B.SubMatrixOf(A, row, col).Set(C))
func (A *FloatMatrix) SetSubMatrix(row, col int, mat *FloatMatrix) error {
    r, c := mat.Size()
    if r+row > A.Rows() || c+col > A.Cols() {
        s := fmt.Sprintf("(%d+%d, %d+%d) > (%d,%d)\n", r, row, c, col, A.Rows(), A.Cols())
        return errors.New(s)
    }
    for i := 0; i < r; i++ {
        for j := 0; j < c; j++ {
            A.SetAt(row+i, col+j, mat.GetAt(i, j))
        }
    }
    return nil
}

// Get sub-matrix starting at (row, col). Sizes parameters define (nrows, ncols) number of
// rows and number of columns. If len(sizes) is zero size is then (nrows, ncols) is
// (Rows()-row, Cols()-col).If len(sizes) is one then (nrows, ncols) is (sizes[0], Cols()-col)
// In all other cases (nrows, ncols) is (sizes[0], sizes[1]). 
// Return nil if nrows+row >= A.Rows() or ncols+col >= A.Cols()
// **DEPRECEATED**
func (A *FloatMatrix) GetSubMatrix(row, col int, sizes ...int) (m *FloatMatrix) {
    var nrows, ncols int = 0, 0
    switch len(sizes) {
    case 0:
        nrows = A.Rows() - row
        ncols = A.Cols() - col
    case 1:
        nrows = sizes[0]
        ncols = A.Cols() - col
    default:
        nrows = sizes[0]
        ncols = sizes[1]
    }
    if row+nrows > A.Rows() || col+ncols > A.Cols() {
        return nil
    }
    var colArray []float64 = nil
    m = FloatZeros(nrows, ncols)
    for i := 0; i < ncols; i++ {
        colArray = A.GetColumnArray(col+i, colArray)
        m.SetColumnArray(i, colArray[row:])
    }
    return m
}

// Compute A = fn(C) by applying function fn to all elements in indexes.
// For all i in indexes: A[i] = fn(C[i]).
// If C is nil then computes inplace A = fn(A). If C is not nil then sizes of A and C must match.
// Returns pointer to self. **DEPRECEATED**
func (A *FloatMatrix) ApplyToIndexes(C *FloatMatrix, indexes []int, fn func(float64) float64) *FloatMatrix {
    if C != nil && !A.SizeMatch(C.Size()) {
        return nil
    }
    B := C
    if C == nil {
        B = A
    }
    if len(indexes) > 0 {
        nrows := A.Rows()
        step := A.LeadingIndex()
        N := A.NumElements()
        for _, k := range indexes {
            k = (k + N) % N
            rk := realIndex(k, nrows, step)
            A.elements[rk] = fn(B.elements[rk])
        }
    }
    return A
}

// Local Variables:
// tab-width: 4
// End:
