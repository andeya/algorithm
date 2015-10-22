// Copyright (c) Harri Rautila, 2012

// This file is part of go.opt/matrix package. It is free software, distributed
// under the terms of GNU Lesser General Public License Version 3, or any later
// version. See the COPYING tile included in this archive.

package matrix

import "math/cmplx"

// Compute in-place product A[i,j] *= alpha
func (A *ComplexMatrix) Scale(alpha complex128, indexes ...int) *ComplexMatrix {
    nrows := A.Rows()
    step := A.LeadingIndex()
    if len(indexes) == 0 {
        for k := 0; k < A.NumElements(); k++ {
            rk := realIndex(k, nrows, step)
            A.elements[rk] *= alpha
        }
    } else {
        N := A.NumElements()
        for k := 0; k < A.NumElements(); k++ {
            k = (k + N) % N
            rk := realIndex(k, nrows, step)
            A.elements[rk] *= alpha
        }
    }
    return A
}

// Compute in-place A[indexes[i]] *= values[i]. Indexes are in column-major order.
func (A *ComplexMatrix) ScaleIndexes(indexes []int, values []complex128) *ComplexMatrix {
    nrows := A.Rows()
    step := A.LeadingIndex()
    if len(indexes) == 0 {
        return A
    }
    N := A.NumElements()
    for i, k := range indexes {
        if i >= len(values) {
            return A
        }
        k = (k + N) % N
        rk := realIndex(k, nrows, step)
        A.elements[rk] *= values[i]
    }
    return A
}

// Compute in-place sum A[i,j] += alpha
func (A *ComplexMatrix) Add(alpha complex128, indexes ...int) *ComplexMatrix {
    nrows := A.Rows()
    step := A.LeadingIndex()
    N := A.NumElements()
    if len(indexes) == 0 {
        for k := 0; k < N; k++ {
            rk := realIndex(k, nrows, step)
            A.elements[rk] += alpha
        }
    } else {
        for _, k := range indexes {
            k = (k + N) % N
            rk := realIndex(k, nrows, step)
            A.elements[rk] += alpha
        }
    }
    return A
}

// Compute in-place A[indexes[i]] += values[i]. Indexes are in column-major order.
func (A *ComplexMatrix) AddIndexes(indexes []int, values []complex128) *ComplexMatrix {
    if len(indexes) == 0 {
        return A
    }
    nrows := A.Rows()
    step := A.LeadingIndex()
    N := A.NumElements()
    for i, k := range indexes {
        if i >= len(values) {
            return A
        }
        k = (k + N) % N
        rk := realIndex(k, nrows, step)
        A.elements[rk] += values[i]
    }
    return A
}

// Compute in-place inverse A[i,j] = 1.0/A[i,j]. If indexes is empty calculates for
// all elements
func (A *ComplexMatrix) Inv(indexes ...int) *ComplexMatrix {
    nrows := A.Rows()
    step := A.LeadingIndex()
    N := A.NumElements()
    if len(indexes) == 0 {
        for k := 0; k < N; k++ {
            rk := realIndex(k, nrows, step)
            A.elements[rk] = 1.0 / A.elements[rk]
        }
    } else {
        for _, k := range indexes {
            k = (k + N) % N
            rk := realIndex(k, nrows, step)
            A.elements[rk] = 1.0 / A.elements[rk]
        }
    }
    return A
}

// Compute in-place negation -A[i,j]
func (A *ComplexMatrix) Neg() *ComplexMatrix {
    for k, v := range A.elements {
        A.elements[k] = -v
    }
    return A
}

// Compute element-wise division A /= B. Return A. If A and B sizes
// do not match A is returned unaltered.
func (A *ComplexMatrix) Div(B *ComplexMatrix) *ComplexMatrix {
    if !A.SizeMatch(B.Size()) {
        return A
    }
    N := A.NumElements()
    for k := 0; k < N; k++ {
        rka := realIndex(k, A.Rows(), A.LeadingIndex())
        rkb := realIndex(k, B.Rows(), B.LeadingIndex())
        A.elements[rka] /= B.elements[rkb]
    }
    return A
}

// Compute element-wise product A *= B. Return A. If A and B sizes
// do not match A is returned unaltered.
func (A *ComplexMatrix) Mul(B *ComplexMatrix) *ComplexMatrix {
    if !A.SizeMatch(B.Size()) {
        return A
    }
    N := A.NumElements()
    for k := 0; k < N; k++ {
        rka := realIndex(k, A.Rows(), A.LeadingIndex())
        rkb := realIndex(k, B.Rows(), B.LeadingIndex())
        A.elements[rka] *= B.elements[rkb]
    }
    return A
}

// Compute element-wise sum A += B. Return A. If A and B sizes
// do not match A is returned unaltered.
func (A *ComplexMatrix) Plus(B *ComplexMatrix) *ComplexMatrix {
    if !A.SizeMatch(B.Size()) {
        return A
    }
    N := A.NumElements()
    for k := 0; k < N; k++ {
        rka := realIndex(k, A.Rows(), A.LeadingIndex())
        rkb := realIndex(k, B.Rows(), B.LeadingIndex())
        A.elements[rka] += B.elements[rkb]
    }
    return A
}

// Compute element-wise difference A -= B. Return A. If A and B sizes
// do not match A is returned unaltered.
func (A *ComplexMatrix) Minus(B *ComplexMatrix) *ComplexMatrix {
    if !A.SizeMatch(B.Size()) {
        return A
    }
    N := A.NumElements()
    for k := 0; k < N; k++ {
        rka := realIndex(k, A.Rows(), A.LeadingIndex())
        rkb := realIndex(k, B.Rows(), B.LeadingIndex())
        A.elements[rka] -= B.elements[rkb]
    }
    return A
}

// Compute matrix product C = A * B where A is m*p and B is p*n.
// Returns a new m*n matrix. 
func (A *ComplexMatrix) Times(B *ComplexMatrix) *ComplexMatrix {
    if A.Cols() != B.Rows() {
        return nil
    }
    rows := A.Rows()
    cols := B.Cols()
    acols := A.Cols()
    C := ComplexZeros(rows, cols)
    Ar := A.elements
    Br := B.elements

    // Basic idea:
    // Loop through each matrix always in memory order ie. down each column.

    // cc: start of column j in matrix C (=j*C.Rows)
    // cr: index to C row i in column j  (=j*C.Rows+i)
    // br: index to B, row i in column k (=k*B.Rows+i)
    // ar: index to A, row
    cc := 0
    br := 0
    for j := 0; j < cols; j++ {
        ar := 0
        for k := 0; k < acols; k++ {
            // move C index to first row in current column
            cr := cc
            // beta is value of B[k,j]
            beta := Br[br]
            // zero in B[k,j] does not increment value in C[:,j]
            if beta != 0.0 {
                // C[:,j] += A[:,k]*B[k,j]
                for i := 0; i < rows; i++ {
                    C.elements[cr] += Ar[ar] * beta
                    // move to next row in memory order
                    cr += 1
                    ar += 1
                }
            } else {
                // we skipped all rows in this column, move to start of next column
                ar += rows
            }
            // move to next row in B, here ar points to start of next column in A
            br += 1
        }
        // forward to start of next column in C
        cc += rows
    }
    return C
}

// Compute A = fn(A) by applying function fn element wise to A.
// If indexes array is non-empty function is applied to elements of A
// indexed by the contents of indexes.
func (A *ComplexMatrix) Apply(fn func(complex128) complex128, indexes ...int) *ComplexMatrix {
    nrows := A.Rows()
    step := A.LeadingIndex()
    N := A.NumElements()
    if len(indexes) == 0 {
        for k := 0; k < N; k++ {
            rk := realIndex(k, nrows, step)
            A.elements[rk] = fn(A.elements[rk])
        }
    } else {
        for _, k := range indexes {
            k = (k + N) % N
            rk := realIndex(k, nrows, step)
            A.elements[rk] = fn(A.elements[rk])
        }
    }
    return A
}

// Compute A = fn(A, x) by applying function fn element wise to A.
// If indexes array is non-empty function is applied to elements of A
// indexed by the contents of indexes.
func (A *ComplexMatrix) ApplyConst(x complex128, fn func(complex128, complex128) complex128, indexes ...int) *ComplexMatrix {
    nrows := A.Rows()
    step := A.LeadingIndex()
    N := A.NumElements()
    if len(indexes) == 0 {
        for k := 0; k < N; k++ {
            rk := realIndex(k, nrows, step)
            A.elements[rk] = fn(A.elements[rk], x)
        }
    } else {
        for _, k := range indexes {
            k = (k + N) % N
            rk := realIndex(k, nrows, step)
            A.elements[rk] = fn(A.elements[rk], x)
        }
    }
    return A
}

// Compute A = fn(A, x) by applying function fn element wise to A.
//  For all i in indexes: A[indexes[i]] = fn(A[indexes[i]], values[i])
func (A *ComplexMatrix) ApplyConstValues(values []complex128, fn func(complex128, complex128) complex128, indexes ...int) *ComplexMatrix {
    nrows := A.Rows()
    step := A.LeadingIndex()
    N := A.NumElements()
    for i, k := range indexes {
        if i > len(values) {
            return A
        }
        k = (k + N) % N
        rk := realIndex(k, nrows, step)
        A.elements[rk] = fn(A.elements[rk], values[i])
    }
    return A
}

// Compute in-place conjugate A[i,j]
func (A *ComplexMatrix) Conj() *ComplexMatrix {
    return A.Apply(cmplx.Conj)
}

// Compute in-place Exp(A)
func (A *ComplexMatrix) Exp() *ComplexMatrix {
    return A.Apply(cmplx.Exp)
}

// Compute in-place Log(A)
func (A *ComplexMatrix) Log() *ComplexMatrix {
    return A.Apply(cmplx.Log)
}

// Compute in-place Log10(A)
func (A *ComplexMatrix) Log10() *ComplexMatrix {
    return A.Apply(cmplx.Log10)
}

// Compute in-place Sqrt(A)
func (A *ComplexMatrix) Sqrt() *ComplexMatrix {
    return A.Apply(cmplx.Sqrt)
}

// Compute in-place Pow(A, x)
func (A *ComplexMatrix) Pow(x complex128) *ComplexMatrix {
    return A.ApplyConst(x, cmplx.Pow)
}

// Local Variables:
// tab-width: 4
// End:
