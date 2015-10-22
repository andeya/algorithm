// Copyright (c) Harri Rautila, 2012

// This file is part of go.opt/matrix package. It is free software, distributed
// under the terms of GNU Lesser General Public License Version 3, or any later
// version. See the COPYING tile included in this archive.

package cmath

import (
	"github.com/henrylee2cn/analysis/matrix"
	"math/cmplx"
)

// Make a copy C of A and compute  C *= alpha for all elements in the matrix if list of indexes
// is empty. Otherwise compute C[i] *= alpha for i in indexes array.
func Scale(A *matrix.ComplexMatrix, alpha complex128, indexes ...int) *matrix.ComplexMatrix {
	C := A.Copy()
	return C.Scale(alpha, indexes...)
}

// Make a copy C of A and compute for all k in indexes: C[k] *= values[k].
// Indexes are in column-major order.  Returns a new matrix
func ScaleAt(A *matrix.ComplexMatrix, values []complex128, indexes []int) *matrix.ComplexMatrix {
	C := A.Copy()
	if len(indexes) == 0 {
		return C
	}
	Cr := C.ComplexArray()
	N := A.NumElements()
	for i, k := range indexes {
		if i >= len(values) {
			return C
		}
		if k < 0 {
			k = N + k
		}
		Cr[k] *= values[i]
	}
	return C
}

// Make a copy C of A and compute inverse C[i] = 1.0/C[i]. If indexes is empty calculates for
// all elements. Returns a new matrix.
func Inv(A *matrix.ComplexMatrix, indexes ...int) *matrix.ComplexMatrix {
	C := A.Copy()
	return C.Inv()
}

// Make a copy C of A and compute C += alpha for all elements in the matrix if list of indexes
// is empty. Otherwise compute C[i] += alpha for indexes in column-major order.
func Add(A *matrix.ComplexMatrix, alpha complex128, indexes ...int) *matrix.ComplexMatrix {
	C := A.Copy()
	return C.Add(alpha, indexes...)
}

// Make copy C of A and compute  C[indexes[i]] +=  values[i]. Indexes are in column-major order.
// Returns a new matrix.
func AddAt(A *matrix.ComplexMatrix, values []complex128, indexes []int) *matrix.ComplexMatrix {
	C := A.Copy()
	if len(indexes) == 0 {
		return C
	}
	Cr := C.ComplexArray()
	N := A.NumElements()
	for i, k := range indexes {
		if i >= len(values) {
			return C
		}
		if k < 0 {
			k = N + k
		}
		Cr[k] += values[i]
	}
	return C
}

// Compute element wise division C[i,j] = A[i,j] / B[i,j]. Returns new matrix.
func Div(A, B *matrix.ComplexMatrix) *matrix.ComplexMatrix {
	if !A.SizeMatch(B.Size()) {
		return nil
	}
	C := A.Copy()
	return C.Div(B)
}

// Compute element-wise product C[i,j] = A[i,j] * B[i,j]. Returns new matrix.
func Mul(A, B *matrix.ComplexMatrix) *matrix.ComplexMatrix {
	if !A.SizeMatch(B.Size()) {
		return nil
	}
	C := A.Copy()
	return C.Mul(B)
}

// Compute element-wise sum C = A + B. Returns a new matrix.
func Plus(matrices ...*matrix.ComplexMatrix) *matrix.ComplexMatrix {
	if len(matrices) <= 1 {
		if len(matrices) == 1 {
			return matrices[0]
		}
		return nil
	}
	A := matrices[0]
	for _, B := range matrices[1:] {
		if !A.SizeMatch(B.Size()) {
			return nil
		}
	}
	C := A.Copy()
	for _, B := range matrices[1:] {
		C.Plus(B)
	}
	return C
}

// Compute element-wise difference C = A - B. Returns a new matrix.
func Minus(matrices ...*matrix.ComplexMatrix) *matrix.ComplexMatrix {
	if len(matrices) <= 1 {
		if len(matrices) == 1 {
			return matrices[0]
		}
		return nil
	}
	A := matrices[0]
	for _, B := range matrices[1:] {
		if !A.SizeMatch(B.Size()) {
			return nil
		}
	}
	C := A.Copy()
	for _, B := range matrices[1:] {
		C.Minus(B)
	}
	return C
}

// Compute matrix product C = A * B where A is m*p and B is p*n.
// Returns a new m*n matrix.
func Times(A, B *matrix.ComplexMatrix) *matrix.ComplexMatrix {
	if A.Cols() != B.Rows() {
		return nil
	}
	return A.Times(B)
}

// Make a copy C of A and apply function fn element wise to C.
// For indexes is not empty then  C[indexes[i]] = fn(C[indexes[i]]).
// Returns a new matrix.
func Apply(A *matrix.ComplexMatrix, fn func(complex128) complex128, indexes ...int) *matrix.ComplexMatrix {
	if A == nil {
		return nil
	}
	C := A.Copy()
	return C.Apply(fn, indexes...)
}

// Make a copy C of A and apply function fn element wise to C.
// For indexes is not empty then  C[indexes[i]] = fn(C[indexes[i]], x).
// Returns a new matrix.
func ApplyConst(A *matrix.ComplexMatrix, x complex128, fn func(complex128, complex128) complex128, indexes ...int) *matrix.ComplexMatrix {
	if A == nil {
		return nil
	}
	C := A.Copy()
	return C.ApplyConst(x, fn, indexes...)
}

// Makes a copy of A and for all elements pointed by the element of the indexes array
// calculates fn(A[k], values[i]) where k is the i'th value in the indexes array.
func ApplyConstValues(A *matrix.ComplexMatrix, values []complex128, fn func(complex128, complex128) complex128, indexes []int) *matrix.ComplexMatrix {
	if A == nil {
		return A
	}
	C := A.Copy()
	return C.ApplyConstValues(values, fn, indexes)
}

// Compute Abs(A), Returns a new float valued matrix.
func Abs(A *matrix.ComplexMatrix) *matrix.FloatMatrix {
	C := matrix.FloatZeros(A.Rows(), A.Cols())
	Cr := C.FloatArray()
	Ar := A.ComplexArray()
	for k, v := range Ar {
		Cr[k] = cmplx.Abs(v)
	}
	return C
}

// Compute element-wise C = Conj(A). Returns a new matrix.
func Conj(A *matrix.ComplexMatrix, indexes ...int) *matrix.ComplexMatrix {
	return Apply(A, cmplx.Conj, indexes...)
}

// Compute element-wise C = Exp(A). Returns a new matrix.
func Exp(A *matrix.ComplexMatrix, indexes ...int) *matrix.ComplexMatrix {
	return Apply(A, cmplx.Exp, indexes...)
}

// Compute element-wise C = Log(A). Returns a new matrix.
func Log(A *matrix.ComplexMatrix, indexes ...int) *matrix.ComplexMatrix {
	return Apply(A, cmplx.Log, indexes...)
}

// Compute element-wise C = Pow(A). Returns a new matrix.
func Pow(A *matrix.ComplexMatrix, exp complex128, indexes ...int) *matrix.ComplexMatrix {
	return ApplyConst(A, exp, cmplx.Pow, indexes...)
}

// Compute element-wise C = Sqrt(A). Returns a new matrix.
func Sqrt(A *matrix.ComplexMatrix, indexes ...int) *matrix.ComplexMatrix {
	return Apply(A, cmplx.Sqrt, indexes...)
}

// Return Real(A).
func Real(A *matrix.ComplexMatrix) *matrix.FloatMatrix {
	C := matrix.FloatZeros(A.Size())
	Ar := A.ComplexArray()
	Cr := C.FloatArray()
	for i, v := range Ar {
		Cr[i] = real(v)
	}
	return C
}

// Return Imag(A).
func Imag(A *matrix.ComplexMatrix) *matrix.FloatMatrix {
	C := matrix.FloatZeros(A.Size())
	Ar := A.ComplexArray()
	Cr := C.FloatArray()
	for i, v := range Ar {
		Cr[i] = imag(v)
	}
	return C
}

// Return Complex(Real, Imag). Return a new matrix.
func Complex(Real, Imag *matrix.FloatMatrix) *matrix.ComplexMatrix {
	if !Real.SizeMatch(Imag.Size()) {
		return nil
	}
	C := matrix.ComplexZeros(Real.Size())
	Rr := Real.FloatArray()
	Ir := Imag.FloatArray()
	Cr := C.ComplexArray()
	for i, _ := range Rr {
		Cr[i] = complex(Rr[i], Ir[i])
	}
	return C
}

// Local Variables:
// tab-width: 4
// End:
