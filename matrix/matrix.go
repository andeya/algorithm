// Copyright (c) Harri Rautila, 2012

// This file is part of go.opt/matrix package. It is free software, distributed
// under the terms of GNU Lesser General Public License Version 3, or any later
// version. See the COPYING tile included in this archive.

// Package matrix implements column major matrices.
package matrix

// Minimal interface for linear algebra packages BLAS/LAPACK
type Matrix interface {
    // The number of rows in this matrix.
    Rows() int
    // The number of columns in this matrix.
    Cols() int
    // The number of elements in this matrix.
    NumElements() int
    // size of the leading index. 
    LeadingIndex() int
    // Matrix in string format.
    String() string
    // Make a copy  and return as Matrix interface type.
    MakeCopy() Matrix
    // Match size. Return true if equal.
    SizeMatch(int, int) bool
    // Get matrix size. Return pair (rows, cols).
    Size() (int, int)
    // Test for type equality.
    EqualTypes(...Matrix) bool
}

// Interface for real and complex scalars.
type Scalar interface {
    Float() float64
    Complex() complex128
}

// Float constant
type FScalar float64

// Return self as float64.
func (self FScalar) Float() float64 { return float64(self) }

// Return complex(self, 0)
func (self FScalar) Complex() complex128 { return complex(float64(self), 0) }

// Complex constant
type CScalar complex128

// Return real(self).
func (self CScalar) Float() float64 { return float64(real(self)) }

// Return self as complex128.
func (self CScalar) Complex() complex128 { return complex128(self) }

// Stacking direction for matrix constructor.
type Stacking int

const StackDown = Stacking(0)
const StackRight = Stacking(1)

// Matrix constructor data order
type DataOrder int

const RowOrder = DataOrder(0)
const ColumnOrder = DataOrder(1)

// Tridiagonal matrix type
type Tridiagonal int

const Symmetric = Tridiagonal(0)
const Lower = Tridiagonal(1)
const Upper = Tridiagonal(2)

// Matrix dimensions, rows, cols and leading index. For column major matrices 
// leading index is equal to row count.
type dimensions struct {
    rows int
    cols int
    // actual offset between leading index
    step int
}

// Return number of rows.
func (A *dimensions) Rows() int {
    if A == nil {
        return 0
    }
    return A.rows
}

// Return number of columns.
func (A *dimensions) Cols() int {
    if A == nil {
        return 0
    }
    return A.cols
}

// Return size of the matrix as rows, cols pair.
func (A *dimensions) Size() (int, int) {
    if A == nil {
        return 0, 0
    }
    return A.rows, A.cols
}

// Set dimensions. Does not affect element allocations. Optional 3rd argument
// sets also the row stride of data.
func (A *dimensions) SetSize(nrows, ncols int, step ...int) {
    A.rows = nrows
    A.cols = ncols
    if len(step) > 0 {
        A.step = step[0]
    }
}

// Return the leading index size. Column major matrices it is row count.
func (A *dimensions) LeadingIndex() int {
    return A.step
}

// Return total number of elements.
func (A *dimensions) NumElements() int {
    if A == nil {
        return 0
    }
    return A.rows * A.cols
}

// Return true if size of A is equal to parameter size (rows, cols).
func (A *dimensions) SizeMatch(rows, cols int) bool {
    return A != nil && A.rows == rows && A.cols == cols
}

// Change matrix shape if number of elements match to rows*cols. Optional 3rd
// arguments is the new row stride of data. Row stride is not changed unless given
// or if current row stride is equal to number of rows in matrix. On the second case
// new row stride is set equal to new row count.
func Reshape(m Matrix, rows, cols int, steps ...int) {
    step := m.LeadingIndex()
    if m.Rows() == step {
        // this not a submatrix, (not definitive?)
        step = rows
    }
    if len(steps) > 0 {
        step = steps[0]
    }
    if rows*cols == m.NumElements() {
        switch m.(type) {
        case *FloatMatrix:
            m.(*FloatMatrix).SetSize(rows, cols, step)
        case *ComplexMatrix:
            m.(*ComplexMatrix).SetSize(rows, cols, step)
        }
    }
}

// Set x = y ie. copy y to x. Matrices must have same number of elements but are
// not required to have same shape. **DEPRECEATED**
func Set(x, y Matrix) {
    if x.NumElements() != y.NumElements() {
        return
    }
    if !x.EqualTypes(y) {
        return
    }
    switch x.(type) {
    case *FloatMatrix:
        x.(*FloatMatrix).Set(y.(*FloatMatrix))
    case *ComplexMatrix:
        x.(*ComplexMatrix).Set(y.(*ComplexMatrix))
    }
}

// Create a set of indexes from start to end-1 with interval step.
func MakeIndexSet(start, end, step int) []int {
    return Indexes(start, end, step)
    /*
    	if start < 0 {
    		start = 0
    	}
    	if end < 0 {
    		end = 0
    	}
    	if end-start == 0 {
    		return make([]int, 0, 1)
    	}
    	if step < 0 {
    		step = 1
    	}
    	//sz := (end-start)/step + 1
    	inds := make([]int, 0)
    	for k := start; k < end; k += step {
    		inds = append(inds, k)
    	}
    	return inds
    */
}

// Create index set to access diagonal entries of matrix of size (n, n).
func MakeDiagonalSet(n int) []int {
    if n <= 0 {
        return []int{}
    }
    return Indexes(0, n*n, n+1)
}

// Create index set for a row in matrix M. 
func RowIndexes(m Matrix, row int) []int {
    nrows, N := m.Size()
    if row < 0 {
        row += nrows
    }
    if row > nrows {
        return []int{}
    }
    iset := make([]int, N)
    for i := 0; i < N; i++ {
        k := i*N + row
        iset[i] = k
    }
    return iset
}

// Create index set for a column in matrix M. 
func ColumnIndexes(m Matrix, col int) []int {
    N, ncols := m.Size()
    if col < 0 {
        col += ncols
    }
    if col > ncols {
        return []int{}
    }
    iset := make([]int, N)
    for i := 0; i < N; i++ {
        k := col*N + i
        iset[i] = k
    }
    return iset
}

// Create index set for diagonal in matrix M. 
func DiagonalIndexes(m Matrix) []int {
    if m.Cols() != m.Rows() {
        return []int{}
    }
    return Indexes(0, m.NumElements(), m.Rows()+1)
}

// Create an index set. Three argument call is (start, end, step) where start < end and step > 0.
// It produces set of indexes from start to end-1 with step interval. Two argument call is 
// interpreted as (start, end, 1) and one argument call as (0, end, 1).
func Indexes(specs ...int) []int {
    start := 0
    end := 0
    step := 1
    if len(specs) == 1 {
        end = specs[0]
    } else if len(specs) == 2 {
        start = specs[0]
        end = specs[1]
    } else if len(specs) > 2 {
        start = specs[0]
        end = specs[1]
        step = specs[2]
    }
    iset := make([]int, 0)
    // must have: start < end && step > 0
    if !(start < end && step > 0) {
        return iset
    }
    for i := start; i < end; i += step {
        iset = append(iset, i)
    }
    return iset
}

// Translate submatrix direct index to underlying matrix index.
// (relative to the start of submatrix.)
func realIndex(index, nrows, nstep int) int {
    if nrows == nstep {
        return index
    }
    col := index / nrows
    row := index - col*nrows

    return col*nstep + row
}

// Local Variables:
// tab-width: 4
// End:
