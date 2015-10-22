// Copyright (c) Harri Rautila, 2012,2013

// This file is part of go.opt/matrix package. It is free software, distributed
// under the terms of GNU Lesser General Public License Version 3, or any later
// version. See the COPYING tile included in this archive.

package matrix

import "math"

// Compute in-place A *= alpha for all elements in the matrix if list of indexes
// is empty. Otherwise compute A[i] *= alpha for indexes in column-major order.
func (A *FloatMatrix) Scale(alpha float64, indexes ...int) *FloatMatrix {
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
func (A *FloatMatrix) ScaleIndexes(indexes []int, values []float64) *FloatMatrix {
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

// Compute in-place remainder A[i,j] %= alpha
func (A *FloatMatrix) Mod(alpha float64, indexes ...int) *FloatMatrix {
    nrows := A.Rows()
    step := A.LeadingIndex()
    N := A.NumElements()
    if len(indexes) == 0 {
        for k := 0; k < N; k++ {
            rk := realIndex(k, nrows, step)
            A.elements[rk] = math.Mod(A.elements[rk], alpha)
        }
    } else {
        for _, k := range indexes {
            k = (k + N) % N
            rk := realIndex(k, nrows, step)
            A.elements[rk] = math.Mod(A.elements[rk], alpha)
        }
    }
    return A
}

// Compute in-place inverse A[i,j] = 1.0/A[i,j]. If indexes is empty calculates for
// all elements
func (A *FloatMatrix) Inv(indexes ...int) *FloatMatrix {
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

// Compute in-place A += alpha for all elements in the matrix if list of indexes
// is empty. Otherwise compute A[i] += alpha for indexes in column-major order.
func (A *FloatMatrix) Add(alpha float64, indexes ...int) *FloatMatrix {
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
func (A *FloatMatrix) AddIndexes(indexes []int, values []float64) *FloatMatrix {
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

// Compute element-wise division A /= B. Return A. If A and B sizes
// do not match A is returned unaltered.
func (A *FloatMatrix) Div(B *FloatMatrix) *FloatMatrix {
    if !A.SizeMatch(B.Size()) {
        onError("div: size mismatch")
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
func (A *FloatMatrix) Mul(B *FloatMatrix) *FloatMatrix {
    if !A.SizeMatch(B.Size()) {
        onError("mul: size mismatch")
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
func (A *FloatMatrix) Plus(B *FloatMatrix) *FloatMatrix {
    if !A.SizeMatch(B.Size()) {
        onError("plus: size mismatch")
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
// do not match A is returned unaltered. Not guaranteed to work for overlapping
// submatrices.
func (A *FloatMatrix) Minus(B *FloatMatrix) *FloatMatrix {
    if !A.SizeMatch(B.Size()) {
        onError("minus: size mismatch")
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
func (A *FloatMatrix) Times(B *FloatMatrix) *FloatMatrix {
    return Times(A, B)
}

// Compute A = fn(A) by applying function fn element wise to A.
// If indexes array is non-empty function is applied to elements of A
// indexed by the contents of indexes.
func (A *FloatMatrix) Apply(fn func(float64) float64, indexes ...int) *FloatMatrix {
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
//  For all i in A:       A[i]          = fn(A[i], x)           if len(indexes) == 0
//  For all i in indexes: A[indexes[i]] = fn(A[indexes[i]], x)  if len(indexes) > 0
func (A *FloatMatrix) ApplyConst(x float64, fn func(float64, float64) float64, indexes ...int) *FloatMatrix {
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
func (A *FloatMatrix) ApplyConstValues(values []float64, fn func(float64, float64) float64, indexes ...int) *FloatMatrix {
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

// Find element-wise maximum. 
func (A *FloatMatrix) Max(indexes ...int) float64 {
    m := math.Inf(-1)
    nrows := A.Rows()
    step := A.LeadingIndex()
    N := A.NumElements()
    if len(indexes) == 0 {
        for k := 0; k < N; k++ {
            rk := realIndex(k, nrows, step)
            m = math.Max(m, A.elements[rk])
        }
    } else {
        for _, k := range indexes {
            k = (k + N) % N
            rk := realIndex(k, nrows, step)
            m = math.Max(m, A.elements[rk])
        }
    }
    return m
}

// Find element-wise minimum. 
func (A *FloatMatrix) Min(indexes ...int) float64 {
    m := math.Inf(+1)
    nrows := A.Rows()
    step := A.LeadingIndex()
    N := A.NumElements()
    if len(indexes) == 0 {
        for k := 0; k < N; k++ {
            rk := realIndex(k, nrows, step)
            m = math.Min(m, A.elements[rk])
        }
    } else {
        for _, k := range indexes {
            k = (k + N) % N
            rk := realIndex(k, nrows, step)
            m = math.Min(m, A.elements[rk])
        }
    }
    return m
}

// Return sum of elements
func (A *FloatMatrix) Sum(indexes ...int) float64 {
    m := 0.0
    nrows := A.Rows()
    step := A.LeadingIndex()
    N := A.NumElements()
    if len(indexes) == 0 {
        for k := 0; k < N; k++ {
            rk := realIndex(k, nrows, step)
            m += A.elements[rk]
        }
    } else {
        N := A.NumElements()
        for _, k := range indexes {
            k = (k + N) % N
            rk := realIndex(k, nrows, step)
            m += A.elements[rk]
        }
    }
    return m
}

// Compute element-wise A = Exp(A).
func (A *FloatMatrix) Exp() *FloatMatrix {
    return A.Apply(math.Exp)
}

// Compute element-wise A = Sqrt(A).
func (A *FloatMatrix) Sqrt() *FloatMatrix {
    return A.Apply(math.Sqrt)
}

// Compute element-wise A = Log(A).
func (A *FloatMatrix) Log() *FloatMatrix {
    return A.Apply(math.Log)
}

// Compute element-wise A = Pow(A).
func (A *FloatMatrix) Pow(exp float64) *FloatMatrix {
    return A.ApplyConst(exp, math.Pow)
}

// Local Variables:
// tab-width: 4
// End:
