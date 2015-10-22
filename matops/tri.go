// Copyright (c) Harri Rautila, 2013

// This file is part of github.com/hrautila/matops package. It is free software,
// distributed under the terms of GNU Lesser General Public License Version 3, or
// any later version. See the COPYING tile included in this archive.

package matops

import (
	"github.com/henrylee2cn/algorithm/matrix"
)

// Make A tridiagonal, upper, non-unit matrix by clearing the strictly lower part
// of the matrix.
func TriU(A *matrix.FloatMatrix) *matrix.FloatMatrix {
	var Ac matrix.FloatMatrix
	var k int
	mlen := imin(A.Rows(), A.Cols())
	for k = 0; k < mlen; k++ {
		Ac.SubMatrixOf(A, k+1, k, A.Rows()-k-1, 1)
		Ac.SetIndexes(0.0)
	}
	if A.Cols() < A.Rows() {
		Ac.SubMatrixOf(A, A.Cols(), 0)
		Ac.SetIndexes(0.0)
	}
	return A
}

// Make A tridiagonal, upper, unit matrix by clearing the strictly lower part
// of the matrix and setting diagonal elements to one.
func TriUU(A *matrix.FloatMatrix) *matrix.FloatMatrix {
	var Ac matrix.FloatMatrix
	var k int
	mlen := imin(A.Rows(), A.Cols())
	for k = 0; k < mlen; k++ {
		Ac.SubMatrixOf(A, k+1, k, A.Rows()-k-1, 1)
		Ac.SetIndexes(0.0)
		A.SetAt(k, k, 1.0)
	}
	// last element on diagonal
	A.SetAt(k, k, 1.0)
	if A.Cols() < A.Rows() {
		Ac.SubMatrixOf(A, A.Cols(), 0)
		Ac.SetIndexes(0.0)
	}
	return A
}

// Make A tridiagonal, lower, unit matrix by clearing the strictly upper part
// of the matrix and setting diagonal elements to one.
func TriLU(A *matrix.FloatMatrix) *matrix.FloatMatrix {
	var Ac matrix.FloatMatrix
	mlen := imin(A.Rows(), A.Cols())
	A.SetAt(0, 0, 1.0)
	for k := 1; k < mlen; k++ {
		A.SetAt(k, k, 1.0)
		Ac.SubMatrixOf(A, 0, k, k, 1)
		Ac.SetIndexes(0.0)
	}
	if A.Cols() > A.Rows() {
		Ac.SubMatrixOf(A, 0, A.Rows())
		Ac.SetIndexes(0.0)
	}
	return A
}

// Make A tridiagonal, lower, non-unit matrix by clearing the strictly upper part
// of the matrix.
func TriL(A *matrix.FloatMatrix) *matrix.FloatMatrix {
	var Ac matrix.FloatMatrix
	mlen := imin(A.Rows(), A.Cols())
	for k := 1; k < mlen; k++ {
		Ac.SubMatrixOf(A, 0, k, k, 1)
		Ac.SetIndexes(0.0)
	}
	if A.Cols() > A.Rows() {
		Ac.SubMatrixOf(A, 0, A.Rows())
		Ac.SetIndexes(0.0)
	}
	return A
}

// Local Variables:
// tab-width: 4
// indent-tabs-mode: nil
// End:
