// Copyright (c) Harri Rautila, 2012

// This file is part of go.opt/matrix package. It is free software, distributed
// under the terms of GNU Lesser General Public License Version 3, or any later
// version. See the COPYING tile included in this archive.

package matrix

import (
    //"errors"
    //"fmt"
    "math"
    "math/cmplx"
    "math/rand"
)

// A column-major matrix backed by a flat array of all elements.
type FloatMatrix struct {
    dimensions
    // flattened matrix data. elements[i*step+j] is col i, row j
    elements []float64
}

// Interface for producing float numbers
type FloatGenerator interface {
	Next() float64
}

// Produce constant value
type ConstFloat struct {
	Value float64
}

func (g *ConstFloat) Next() float64 {
	return g.Value
}

// Produce uniformly distributed value
type UniformFloat struct {
	Low  float64
	High float64
}

func (g *UniformFloat) Next() float64 {
	return g.Low + (g.High - g.Low) * rand.Float64()
}

// Produce normally distributed value
type NormalFloat struct {
	Mean float64
	StdDev float64
}

func (g *NormalFloat) Next() float64 {
	return g.Mean + g.StdDev * rand.NormFloat64()
}

// Create a column-major matrix from a flat array of elements.
// Assumes values are in column-major order.
func FloatNew(rows, cols int, elements []float64, order ...DataOrder) *FloatMatrix {
    if len(order) > 0 && order[0] != ColumnOrder {
        // data in row order
        table := make([][]float64, 0)
        for i := 0; i < rows; i++ {
            table = append(table, elements[i*cols:(i+1)*cols])
        }
        return FloatMatrixFromTable(table, RowOrder)
    }
    e := make([]float64, rows*cols)
    copy(e, elements)
    return makeFloatMatrix(rows, cols, e)
}

// Create a column major vector from an array of elements. Shorthand for
// call FloatNew(len(elems), 1, elems).
func FloatVector(elements []float64) *FloatMatrix {
    rows := len(elements)
    e := make([]float64, rows)
    copy(e, elements)
    return makeFloatMatrix(rows, 1, e)
}

// Create a singleton matrix from float value. Shorthand for calling
// MakeMatrix(1, 1, value-array-of-length-one).
func FloatValue(value float64) *FloatMatrix {
    e := make([]float64, 1)
    e[0] = value
    return makeFloatMatrix(1, 1, e)
}

// Create random matrix with elements from [0.0, 1.0) uniformly distributed..
func FloatUniform(rows, cols int) *FloatMatrix {
    A := FloatZeros(rows, cols)
    for i, _ := range A.elements {
        A.elements[i] = rand.Float64()
    }
    return A
}

// Create random matrix with elements from normal distribution (mean=0.0, stddev=1.0)
// *DEPRECEATED*
func FloatNormal(rows, cols int) *FloatMatrix {
    A := FloatZeros(rows, cols)
    for i, _ := range A.elements {
        A.elements[i] = rand.NormFloat64()
    }
    return A
}

// Create symmetric n by n random  matrix with elements from [0.0, 1.0).
// *DEPRECEATED*
func FloatUniformSymmetric(n int, uplo ...Tridiagonal) *FloatMatrix {
    var symm Tridiagonal = Symmetric
    if len(uplo) > 0 {
        symm = uplo[0]
    }
    A := FloatZeros(n, n)
    for i := 0; i < n; i++ {
        for j := i; j < n; j++ {
            val := rand.Float64()
            if symm == Symmetric || symm == Upper {
                A.SetAt(i, j, val)
            }
			if symm == Lower {
                A.SetAt(j, i, val)
			}
            if symm == Symmetric && i != j {
                A.SetAt(j, i, val)
            }
        }
    }
    return A
}

// Create symmetric n by n random  matrix with elements from normal distribution.
// *DEPRECEATED*
func FloatNormalSymmetric(n int, uplo ...Tridiagonal) *FloatMatrix {
    var symm Tridiagonal = Symmetric
    if len(uplo) > 0 {
        symm = uplo[0]
    }
    A := FloatZeros(n, n)
    for i := 0; i < n; i++ {
        for j := i; j < n; j++ {
            val := rand.NormFloat64()
            if symm == Symmetric || symm == Upper {
                A.SetAt(i, j, val)
            }
			if symm == Lower {
                A.SetAt(j, i, val)
			}
            if symm == Symmetric && i != j {
                A.SetAt(j, i, val)
            }
        }
    }
    return A
}

// Create a column-major matrix from a array of arrays. Parameter order
// indicates if data is array of rows (RowOrder) or array of columns (ColumnOrder).
func FloatMatrixFromTable(data [][]float64, order ...DataOrder) *FloatMatrix {
    var rows, cols int
    if len(order) == 0 || order[0] == ColumnOrder {
        cols = len(data)
        rows = len(data[0])
    } else {
        rows = len(data)
        cols = len(data[0])
    }

    if rows*cols == 0 {
        return FloatZeros(rows, cols)
    }
    elements := make([]float64, rows*cols)
    if len(order) == 0 || order[0] == ColumnOrder {
        for i := 0; i < cols; i++ {
            copy(elements[i*rows:], data[i][0:])
        }
    } else {
        for i := 0; i < cols; i++ {
            for j := 0; j < rows; j++ {
                elements[i*rows+j] = data[j][i]
            }
        }
    }
    return makeFloatMatrix(rows, cols, elements)
}

// Create a new matrix from a list of matrices. New matrix has dimension (M, colmax)
// if direction is StackDown, and (rowmax, N) if direction is StackRight. 
// M is sum of row counts of argument matrices and N is sum of column counts of arguments.
// Colmax is the largest column count of matrices and rowmax is the largest row count.
// Return  new matrix and array of submatrix sizes, row counts for StackDown and column
// counts for StackRight
func FloatMatrixStacked(direction Stacking, mlist ...*FloatMatrix) (*FloatMatrix, []int) {
    maxc := 0
    maxr := 0
    N := 0
    M := 0
    for _, m := range mlist {
        m, n := m.Size()
        M += m
        N += n
        if m > maxr {
            maxr = m
        }
        if n > maxc {
            maxc = n
        }
    }
    var mat *FloatMatrix
    indexes := make([]int, 0)
    if direction == StackDown {
        mat = FloatZeros(M, maxc)
        row := 0
        for _, m := range mlist {
            mat.SetSubMatrix(row, 0, m)
            indexes = append(indexes, m.Rows())
            row += m.Rows()
        }
    } else {
        mat = FloatZeros(maxr, N)
        col := 0
        for _, m := range mlist {
            mat.SetSubMatrix(0, col, m)
            indexes = append(indexes, m.Cols())
            col += m.Cols()
        }
    }
    return mat, indexes
}

// Create new zero filled matrix.
func FloatZeros(rows, cols int) *FloatMatrix {
    A := makeFloatMatrix(rows, cols, make([]float64, rows*cols))
    return A
}

// Create new matrix initialized to one.
// *DEPRECEATED*
func FloatOnes(rows, cols int) *FloatMatrix {
    return FloatWithValue(rows, cols, 1.0)
}

// Create new matrix initialized to value.
func FloatWithValue(rows, cols int, value float64) *FloatMatrix {
    A := FloatZeros(rows, cols)
    for k, _ := range A.elements {
        A.elements[k] = value
    }
    return A
}

// Create new identity matrix. Row count must equal column count.
// *DEPRECEATED*
func FloatIdentity(rows int) *FloatMatrix {
    return FloatDiagonal(rows, 1.0)
}

// Make a square matrix with diagonal set to values. If len(values) is one
// then all entries on diagonal is set to values[0]. If len(values) is
// greater than one then diagonals are set from the list values starting
// from (0,0) until the diagonal is full or values are exhausted.
// *DEPRECEATED* 
func FloatDiagonal(rows int, values ...float64) *FloatMatrix {
    A := FloatZeros(rows, rows)
    step := A.LeadingIndex()
    if len(values) == 1 {
        for k := 0; k < rows; k++ {
            A.elements[k*step+k] = values[0]
        }
    } else {
        for k := 0; k < rows && k < len(values); k++ {
            A.elements[k*step+k] = values[k]
        }
    }
    return A
}

// Make B a submatrix of A with top left corner at (row, col).
//
// If no size argument is given then size of (A.Rows()-row, A.Cols()-col) is assumed.
// The variadic size argument can be either (nrows, ncols) or (nrows, ncols, nstep). 
// The first form sets B size to (nrows, ncols). The second form is used to create
// diagonal vectors over A by setting nrows to one, ncols properly and nstep to
// A.LeadingIndex()+1. Other combinations of parameters in three argument form may
// create unexpected access patterns to underlying matrix.
func (A *FloatMatrix) SubMatrix(B *FloatMatrix, row, col int, size ...int) *FloatMatrix {
    nrows := A.Rows() - row
    ncols := A.Cols() - col
    nstep := A.LeadingIndex()
    switch {
    case len(size) == 2:
        nrows = size[0]
        ncols = size[1]
    case len(size) > 2:
        nrows = size[0]
        ncols = size[1]
        nstep = size[2]
    }
	offset := col*A.LeadingIndex()+row
	if offset < len(A.elements) {
		B.elements = A.elements[offset:]
	} else {
		B.elements = nil
	}
    B.rows = nrows
    B.cols = ncols
    B.step = nstep
	return B
}

// Make B diagonal of matrix A as submatrix vector.
func (A *FloatMatrix) Diag(B *FloatMatrix) *FloatMatrix {
    return A.SubMatrix(B, 0, 0, 1, A.Rows(), A.LeadingIndex()+1)
}


// Set A to be submatrix of B starting from position row, col. Returns A.
// The size argument can be either: nrows, ncols or nrows, ncols, nstep.
// The first form is used to create row or column vectors or normal submatrices.
// Three argument form is used to create diagonal vectors by setting nrows to one,
// ncols to number of columns and nstep to A.LeadingIndex()+1. Other combinations
// of parameters in three argument form may create unexpected access patterns to
// underlying matrix. 
func (A *FloatMatrix) SubMatrixOf(B *FloatMatrix, row, col int, size ...int) *FloatMatrix {
    nrows := B.Rows() - row
    ncols := B.Cols() - col
    nstep := B.LeadingIndex()
    switch {
    case len(size) == 2:
        nrows = size[0]
        ncols = size[1]
    case len(size) > 2:
        nrows = size[0]
        ncols = size[1]
        nstep = size[2]
    }
	offset := col*B.LeadingIndex()+row
	if offset < len(B.elements) {
		A.elements = B.elements[offset:]
	} else {
		A.elements = nil
	}
    A.rows = nrows
    A.cols = ncols
    A.step = nstep
    return A
}

// Make A diagonal of matrix B as submatrix vector.
func (A *FloatMatrix) DiagOf(B *FloatMatrix) *FloatMatrix {
    if B.Rows() < B.Cols() {
        return A.SubMatrixOf(B, 0, 0, 1, B.Rows(), B.LeadingIndex()+1)
    } else if B.Rows() > B.Cols() {
        return A.SubMatrixOf(B, 0, 0, 1, B.Cols(), B.LeadingIndex()+1)
    }
    // here B.Rows() == B.Cols(): standard square matrix
    return A.SubMatrixOf(B, 0, 0, 1, B.Rows(), B.LeadingIndex()+1)
}

// Return the flat column-major element array.
func (A *FloatMatrix) FloatArray() []float64 {
    if A == nil {
        return nil
    }
    return A.elements
}

// Return nil for complex array 
func (A *FloatMatrix) ComplexArray() []complex128 {
    return nil
}

// Return the first element column-major element array.
func (A *FloatMatrix) Float() float64 {
    if A == nil {
        return math.NaN()
    }
    return A.elements[0]
}

// Return Nan for complex singleton.
func (A *FloatMatrix) Complex() complex128 {
    return cmplx.NaN()
}

// Test if parameter matrices are of same type as self.
func (A *FloatMatrix) EqualTypes(mats ...Matrix) bool {
loop:
    for _, m := range mats {
        if m == nil {
            continue loop
        }
        switch m.(type) {
        case *FloatMatrix: // of same type, NoOp
        default: // all others fail.
            return false
        }
    }
    return true
}

// Get the element in the i'th row and j'th column.
func (A *FloatMatrix) GetAt(i int, j int) (val float64) {
    step := A.LeadingIndex()
    //val = A.elements[j*step:j*step+A.Cols()][i]
    if i < 0 {
        i += A.Rows()
    }
    if j < 0 {
        j += A.Cols()
    }
    val = A.elements[j*step+i]
    return
}

// Get elements from column-major indexes. Return new array.
func (A *FloatMatrix) GetIndexes(indexes ...int) []float64 {
    vals := make([]float64, 0)
    N := A.NumElements()
    for _, k := range indexes {
        k = (k + N) % N
        rk := realIndex(k, A.Rows(), A.LeadingIndex())
        vals = append(vals, A.elements[rk])
    }
    return vals
}

// Get element value from column-major index position.
func (A *FloatMatrix) GetIndex(i int) float64 {
    i = (i + A.NumElements()) % A.NumElements()
    rk := realIndex(i, A.Rows(), A.LeadingIndex())
    return A.elements[rk]
}

// Copy A to B, A and B number of elements need not match.
// Copies min(A.NumElements(), B.NumElements()) from start of A to start of B.
func (A *FloatMatrix) CopyTo(B *FloatMatrix) error {
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

// Set A = B, copy values, A and B sizes must match.
func (A *FloatMatrix) Set(B *FloatMatrix) error {
    if !A.SizeMatch(B.Size()) {
        return onError("A != B: size mismatch")
    }
    N := A.NumElements()
    for k := 0; k < N; k++ {
        rka := realIndex(k, A.Rows(), A.LeadingIndex())
        rkb := realIndex(k, B.Rows(), B.LeadingIndex())
        A.elements[rka] = B.elements[rkb]
    }
    return nil
}

// Set the element in the i'th row and j'th column to val.
func (A *FloatMatrix) SetAt(i, j int, val float64) {
    step := A.LeadingIndex()
    if i < 0 {
        i = A.Rows() + i
    }
    if j < 0 {
        j = A.Cols() + j
    }
    A.elements[j*step+i] = val
}

// Set value of singleton matrix.
func (A *FloatMatrix) SetValue(val float64) {
    A.elements[0] = val
}

// Set element values in column-major ordering. Negative indexes are relative 
// to the last element of the matrix. If len(indexes) is zero sets all elements.
func (A *FloatMatrix) SetIndexes(val float64, indexes ...int) {
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

// Set value of i'th element. 
func (A *FloatMatrix) SetIndex(i int, val float64) {
    A.SetIndexes(val, i)
}

// Set values of indexed elements. 
func (A *FloatMatrix) SetIndexesFromArray(values []float64, indexes ...int) {
    nrows := A.Rows()
    nstep := A.LeadingIndex()
    N := A.NumElements()
    for i, k := range indexes {
        if i >= len(values) {
            break
        }
        k = (k + N) % N
        rk := realIndex(k, nrows, nstep)
        A.elements[rk] = values[i]
    }
}

func (A *FloatMatrix) SetFrom(g FloatGenerator, indexes ...int) *FloatMatrix {
	if len(indexes) == 0 {
		for k := 0; k < A.NumElements(); k++ {
			rk := realIndex(k, A.Rows(), A.LeadingIndex())
			A.elements[rk] = g.Next()
		}
	} else {
		for _, k := range indexes {
			rk := realIndex(k, A.Rows(), A.LeadingIndex())
			A.elements[rk] = g.Next()
		}
	}
	return A
}

func (A *FloatMatrix) SetFromTrm(g FloatGenerator, tridiagonal ...Tridiagonal) *FloatMatrix {
	which := Symmetric
	if len(tridiagonal) > 0 {
		which = tridiagonal[0]
	}
	for i := 0; i < A.Rows(); i++ {
		A.SetAt(i, i, g.Next())
		switch (which) {
		case Lower:
			for j := 0; j < i; j++ {
				A.SetAt(i, j, g.Next())
			}
		case Upper:
			for j := i+1; j < A.Cols(); j++ {
				A.SetAt(i, j, g.Next())
			}
		default: // Symmetric
			for j := i+1; j < A.Cols(); j++ {
				v := g.Next()
				A.SetAt(i, j, v)
				A.SetAt(j, i, v)
			}
		}
	}
	return A
}

// Create a copy of matrix.
func (A *FloatMatrix) Copy() (B *FloatMatrix) {
    B = new(FloatMatrix)
    B.elements = make([]float64, A.NumElements())
    B.SetSize(A.Rows(), A.Cols(), A.Rows())
    B.Set(A)
    return
}

func (A *FloatMatrix) MakeCopy() Matrix {
    return A.Copy()
}

// Copy and transpose matrix. Returns new matrix.
func (A *FloatMatrix) Transpose() *FloatMatrix {
    rows := A.Rows()
    cols := A.Cols()
    newelems := transposeFloatArray(rows, cols, A.LeadingIndex(), A.elements)
    return makeFloatMatrix(cols, rows, newelems)
}

// Transpose a column major data array.
func transposeFloatArray(rows, cols, step int, data []float64) []float64 {
    newelems := make([]float64, rows*cols)
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
func makeFloatMatrix(rows, cols int, elements []float64) *FloatMatrix {
    A := new(FloatMatrix)
    A.SetSize(rows, cols, rows)
    A.elements = elements
    return A
}


// Local Variables:
// tab-width: 4
// End:
