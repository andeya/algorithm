// Copyright (c) Harri Rautila, 2012

// This file is part of github.com/hrautila/matrix package.
// It is free software, distributed under the terms of
// GNU Lesser General Public License Version 3, or any later
// version. See the COPYING tile included in this archive.

package matrix

import (
    "errors"
    "math"
    "math/rand"
)

// A column-major dense matrix backed by a flat array of all elements.
type ComplexMatrix struct {
    dimensions
    // flattened matrix data. elements[i*step+j] is col i, row j
    elements []complex128
}

// Create a column-major matrix from a flat array of elements.
func ComplexNew(rows, cols int, elements []complex128) *ComplexMatrix {
    e := make([]complex128, rows*cols)
    copy(e, elements)
    return makeComplexMatrix(rows, cols, e)
}

// Create a column major vector from an array of elements
func ComplexVector(elements []complex128) *ComplexMatrix {
    rows := len(elements)
    e := make([]complex128, rows)
    copy(e, elements)
    return makeComplexMatrix(rows, 1, e)
}

// Create a singleton matrix from flot value.
func ComplexValue(value complex128) *ComplexMatrix {
    e := make([]complex128, 1)
    e[0] = value
    return makeComplexMatrix(1, 1, e)
}

// Create random matrix with element's real and imaginary parts
// from [0.0, 1.0).
func ComplexUniform(rows, cols int) *ComplexMatrix {
    A := ComplexZeros(rows, cols)
    for i, _ := range A.elements {
        re := rand.Float64()
        im := rand.Float64()
        A.elements[i] = complex(re, im)
    }
    return A
}

// Create symmetric n by n random  matrix with element's real and imaginary
// parts from [0.0, 1.0).
func ComplexUniformSymmetric(n int) *ComplexMatrix {
    A := ComplexZeros(n, n)
    for i := 0; i < n; i++ {
        for j := i; j < n; j++ {
            re := rand.Float64()
            im := rand.Float64()
            val := complex(re, im)
            A.SetAt(val, i, j)
            if i != j {
                A.SetAt(val, j, i)
            }
        }
    }
    return A
}

// Create random matrix with element's real and imaginary parts
// from [0.0, 1.0).
func ComplexNormal(rows, cols int) *ComplexMatrix {
    A := ComplexZeros(rows, cols)
    for i, _ := range A.elements {
        re := rand.NormFloat64()
        im := rand.NormFloat64()
        A.elements[i] = complex(re, im)
    }
    return A
}

// Create symmetric n by n random  matrix with element's real and imaginary
// parts from normal distribution.
func ComplexNormalSymmetric(n int) *ComplexMatrix {
    A := ComplexZeros(n, n)
    for i := 0; i < n; i++ {
        for j := i; j < n; j++ {
            re := rand.NormFloat64()
            im := rand.NormFloat64()
            val := complex(re, im)
            A.SetAt(val, i, j)
            if i != j {
                A.SetAt(val, j, i)
            }
        }
    }
    return A
}

// Create a column-major matrix from a array of arrays. Parameter rowOrder
// indicates if data is array of rows or array of columns.
func ComplexMatrixFromTable(data [][]complex128, order DataOrder) *ComplexMatrix {
    var rows, cols int
    if order == RowOrder {
        rows = len(data)
        cols = len(data[0])
    } else {
        cols = len(data)
        rows = len(data[0])
    }
    elements := make([]complex128, rows*cols)
    if order == RowOrder {
        for i := 0; i < cols; i++ {
            for j := 0; j < rows; j++ {
                elements[i*rows+j] = data[j][i]
            }
        }
    } else {
        for i := 0; i < cols; i++ {
            copy(elements[i*rows:], data[i][0:])
        }
    }
    return makeComplexMatrix(rows, cols, elements)
}

// Create new zero filled matrix.
func ComplexZeros(rows, cols int) *ComplexMatrix {
    A := makeComplexMatrix(rows, cols, make([]complex128, rows*cols))
    return A
}

// Create new matrix initialized to one.
func ComplexOnes(rows, cols int) *ComplexMatrix {
    return ComplexNumbers(rows, cols, complex(1.0, 0.0))
}

// Create new matrix initialized to value.
func ComplexNumbers(rows, cols int, value complex128) *ComplexMatrix {
    A := ComplexZeros(rows, cols)
    for k, _ := range A.elements {
        A.elements[k] = value
    }
    return A
}

// Create new identity matrix. Row count must equal column count.
func ComplexIdentity(rows, cols int) (A *ComplexMatrix, err error) {
    A = nil
    if rows != cols {
        err = ErrorDimensionMismatch
        return
    }
    A = ComplexZeros(rows, cols)
    step := A.LeadingIndex()
    for k := rows; k < rows; k++ {
        A.elements[k*step+k] = complex(1.0, 0.0)
    }
    return
}

// Make a submatrix of A starting from position row, col. Returns a new matrix.
// If size not given it is assumed to be [A.Rows()-row, A.Cols()-col].
func (A *ComplexMatrix) SubMatrix(row, col int, size ...int) *ComplexMatrix {
    var nrows, ncols int
    M := new(ComplexMatrix)
    if len(size) < 2 {
        nrows = A.Rows() - row
        ncols = A.Cols() - col
    } else {
        nrows = size[0]
        ncols = size[1]
    }
    // can we support mapping a matrix on a vector ??
    M.elements = A.elements[col*A.LeadingIndex()+row:]
    M.rows = nrows
    M.cols = ncols
    M.step = A.LeadingIndex()
    return M
}

// Set A to be submatrix of B
func (A *ComplexMatrix) SubMatrixOf(B *ComplexMatrix, row, col int, size ...int) *ComplexMatrix {
    nrows := B.Rows() - row
    ncols := B.Cols() - col
    if len(size) >= 2 {
        nrows = size[0]
        ncols = size[1]
    }
    A.elements = B.elements[col*B.LeadingIndex()+row:]
    A.rows = nrows
    A.cols = ncols
    A.step = B.LeadingIndex()
    return A
}

// Return nil for float array 
func (A *ComplexMatrix) FloatArray() []float64 {
    return nil
}

// Return the flat column-major element array.
func (A *ComplexMatrix) ComplexArray() []complex128 {
    return A.elements
}

// Return Nan for float singleton.
func (A *ComplexMatrix) Float() float64 {
    return math.NaN()
}

// Return the first element column-major element array.
func (A *ComplexMatrix) Complex() complex128 {
    return A.elements[0]
}

// Return true for complex matrix.
func (A *ComplexMatrix) IsComplex() bool {
    return true
}

// Test if parameter matrices are of same type as self.
func (A *ComplexMatrix) EqualTypes(mats ...Matrix) bool {
loop:
    for _, m := range mats {
        if m == nil {
            continue loop
        }
        switch m.(type) {
        case *ComplexMatrix: // of same type, NoOp
        default: // all others fail.
            return false
        }
    }
    return true
}

// Get the element in the i'th row and j'th column.
func (A *ComplexMatrix) GetAt(i int, j int) (val complex128) {
    step := A.LeadingIndex()
    val = A.elements[j*step : j*step+A.Cols()][i]
    return
}

// Get elements from column-major indexes. Return new array.
func (A *ComplexMatrix) GetIndexes(indexes ...int) []complex128 {
    vals := make([]complex128, 0)
    N := A.NumElements()
    for _, k := range indexes {
        k = (N + k) % N
        rk := realIndex(k, A.Rows(), A.LeadingIndex())
        vals = append(vals, A.elements[rk])
    }
    return vals
}

// Get i'th element in column-major ordering
func (A *ComplexMatrix) GetIndex(i int) complex128 {
    if i < 0 {
        i = A.NumElements() + i
    }
    i %= A.NumElements()
    return A.elements[i]
}

// Get values for indexed elements.  **DEPRECEATED**
func (A *ComplexMatrix) GetIndexesFromArray(indexes []int) []complex128 {
    vals := make([]complex128, len(indexes))
    for i, k := range indexes {
        if k < 0 {
            k = A.NumElements() + k
        }
        k %= A.NumElements()
        vals[i] = A.elements[k]
    }
    return vals
}

// Get copy of i'th row.  **DEPRECEATED**
func (A *ComplexMatrix) GetRowArray(i int, vals []complex128) []complex128 {
    if cap(vals) < A.Cols() {
        vals = make([]complex128, A.Cols())
    }
    step := A.LeadingIndex()
    for j := 0; j < A.Cols(); j++ {
        vals[j] = A.elements[j*step+i]
    }
    return vals
}

// Get copy of i'th column.  **DEPRECEATED**
func (A *ComplexMatrix) GetColumnArray(i int, vals []complex128) []complex128 {
    if cap(vals) < A.Rows() {
        vals = make([]complex128, A.Rows())
    }
    step := A.LeadingIndex()
    for j := 0; j < A.Rows(); j++ {
        vals[j] = A.elements[i*step+j]
    }
    return vals
}

// Get copy of i'th row. Return parameter matrix. If vec is too small 
// reallocate new vector and return it.  **DEPRECEATED**
func (A *ComplexMatrix) GetRow(i int, vec *ComplexMatrix) *ComplexMatrix {
    if vec == nil || vec.NumElements() < A.Cols() {
        vec = ComplexZeros(A.Cols(), 1)
    }
    step := A.LeadingIndex()
    ar := vec.ComplexArray()
    for j := 0; j < A.Cols(); j++ {
        ar[j] = A.elements[j*step+i]
    }
    return vec
}

// Get copy of i'th column. See GetRow.  **DEPRECEATED**
func (A *ComplexMatrix) GetColumn(i int, vec *ComplexMatrix) *ComplexMatrix {
    if vec == nil || vec.NumElements() < A.Rows() {
        vec = ComplexZeros(A.Rows(), 1)
    }
    step := A.LeadingIndex()
    ar := vec.ComplexArray()
    for j := 0; j < A.Rows(); j++ {
        ar[j] = A.elements[i*step+j]
    }
    return vec
}

// Get a slice from the underlying storage array. Changing entries
// in the returned slices changes the matrix. Be carefull with this.
// **DEPRECEATED**
func (A *ComplexMatrix) GetSlice(start, end int) []complex128 {
    if start < 0 {
        start = 0
    }
    if end > A.NumElements() {
        end = A.NumElements()
    }
    return A.elements[start:end]
}

// Set A = B, copy values, A and B sizes must match.
func (A *ComplexMatrix) Set(B *ComplexMatrix) error {
    if !A.SizeMatch(B.Size()) {
        return errors.New("A != B: size mismatch")
    }
    ldB := B.LeadingIndex()
    ldA := A.LeadingIndex()
    nrows := A.Rows()
    for k := 0; k < A.Cols(); k++ {
        copy(A.elements[k*ldA:], B.elements[k*ldB:k*ldB+nrows])
    }
    return nil
}

// Set the element in the i'th row and j'th column to val.
func (A *ComplexMatrix) SetAt(val complex128, i int, j int) {
    step := A.LeadingIndex()
    if i < 0 {
        i = A.Rows() + i
    }
    if j < 0 {
        j = A.Cols() + j
    }
    A.elements[j*step+i] = val
}

// Set element values in column-major ordering. Negative indexes are relative 
// to the last element of the matrix. If len(indexes) is zero sets all elements.
func (A *ComplexMatrix) SetIndexes(val complex128, indexes ...int) {
    nrows := A.Rows()
    nstep := A.LeadingIndex()
    N := A.NumElements()
    if len(indexes) == 0 {
        for k := 0; k < N; k++ {
            rk := realIndex(k, nrows, nstep)
            A.elements[rk] = val
        }
        return
    }
    for _, i := range indexes {
        i = (i + N) % N
        rk := realIndex(i, nrows, nstep)
        A.elements[rk] = val
    }
}

// Set i'th element in column-major ordering
func (A *ComplexMatrix) SetIndex(i int, v complex128) {
    A.SetIndexes(v, i)
}

// Set values of indexed elements.  **DEPRECEATED**
func (A *ComplexMatrix) SetIndexesFromArray(indexes []int, values []complex128) {
    N := A.NumElements()
    for i, k := range indexes {
        if i >= len(values) {
            break
        }
        k = (k + N) % N
        rk := realIndex(k, A.Rows(), A.LeadingIndex())
        A.elements[rk] = values[i]
    }
}

// Set values of i'th row. **DEPRECEATED**
func (A *ComplexMatrix) SetRowArray(i int, vals []complex128) {
    step := A.LeadingIndex()
    for j := 0; j < A.Cols(); j++ {
        A.elements[j*step+i] = vals[j]
    }
}

// Set values of i'th column. **DEPRECEATED**
func (A *ComplexMatrix) SetColumnArray(i int, vals []complex128) {
    step := A.LeadingIndex()
    for j := 0; j < A.Rows(); j++ {
        A.elements[i*step+j] = vals[j]
    }
}

// Create a copy of matrix.
func (A *ComplexMatrix) Copy() (B *ComplexMatrix) {
    B = new(ComplexMatrix)
    B.elements = make([]complex128, A.NumElements())
    B.SetSize(A.Rows(), A.Cols(), A.Rows())
    B.Set(A)
    //copy(B.elements, A.elements)
    return
}

// Create a copy of matrix.
func (A *ComplexMatrix) MakeCopy() Matrix {
    return A.Copy()
}

// Copy A to B, A and B number of elements need not match.
// Copies min(A.NumElements(), B.NumElements()) from start of A to start of B.
func (A *ComplexMatrix) CopyTo(B *ComplexMatrix) error {
    N := A.NumElements()
    if N > B.NumElements() {
        N = B.NumElements()
    }
    for k := 0; k < N; k++ {
        rka := realIndex(k, A.Rows(), A.LeadingIndex())
        rkb := realIndex(k, B.Rows(), B.LeadingIndex())
        B.elements[rkb] = A.elements[rka]
    }
    return nil
}

// Copy and transpose matrix. Returns new matrix.
func (A *ComplexMatrix) Transpose() *ComplexMatrix {
    rows := A.Rows()
    cols := A.Cols()
    newelems := transposeComplexArray(rows, cols, A.LeadingIndex(), A.elements)
    return makeComplexMatrix(cols, rows, newelems)
}

// Transpose a column major data array.
func transposeComplexArray(rows, cols, step int, data []complex128) []complex128 {
    newelems := make([]complex128, rows*cols)
    for i := 0; i < rows; i++ {
        for j := 0; j < cols; j++ {
            curI := j*step + i
            newI := i*cols + j
            //fmt.Printf("r: %d, c: %d, move: %d -> %d\n", i, j, curI, newI)
            newelems[newI] = data[curI]
        }
    }
    return newelems
}

// Create a column-major matrix from a flat array of elements. Elements
// slice is not copied to internal elements but assigned, so underlying
// array holding the actual values stays the same.
func makeComplexMatrix(rows, cols int, elements []complex128) *ComplexMatrix {
    A := new(ComplexMatrix)
    A.SetSize(rows, cols, rows)
    A.elements = elements
    return A
}

/*
func applyTest(A, B Matrix, rfunc func(float64,float64)bool, cfunc func(complex128,complex128)bool) bool {
	return false
}
*/

// Local Variables:
// tab-width: 4
// End:
