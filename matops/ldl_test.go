// Copyright (c) Harri Rautila, 2012,2013

// This file is part of github.com/hrautila/matops package. It is free software,
// distributed under the terms of GNU Lesser General Public License Version 3, or
// any later version. See the COPYING tile included in this archive.

package matops

import (
	"github.com/henrylee2cn/algorithm/matrix"
	"testing"
)

func _TestLDLnoPiv(t *testing.T) {
	N := 42
	nb := 8

	A0 := matrix.FloatUniform(N, N)
	A := matrix.FloatZeros(N, N)
	Mult(A, A0, A0, 1.0, 1.0, TRANSB)

	B := matrix.FloatNormal(A.Rows(), 2)
	w := matrix.FloatWithValue(A.Rows(), 2, 1.0)

	// B0 = A*B
	B0 := B.Copy()

	nb = 2
	L, _ := DecomposeLDLnoPiv(A.Copy(), w, LOWER, nb)
	Mult(B0, A, B, 1.0, 0.0, NOTRANS)
	SolveLDLnoPiv(B0, L, LOWER)
	t.Logf("L*D*L.T: ||B - A*X||_1: %e\n", NormP(B0.Minus(B), NORM_ONE))

	U, _ := DecomposeLDLnoPiv(A.Copy(), w, UPPER, nb)
	Mult(B0, A, B, 1.0, 0.0, NOTRANS)
	SolveLDLnoPiv(B0, U, UPPER)
	t.Logf("U*D*U.T: ||B - A*X||_1: %e\n", NormP(B0.Minus(B), NORM_ONE))

}

func TestLDLlower(t *testing.T) {
	/*
	   Ldata := [][]float64{
	    []float64{7.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	    []float64{7.0, 6.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	    []float64{7.0, 6.0, 5.0, 0.0, 0.0, 0.0, 0.0},
	    []float64{7.0, 6.0, 5.0, 4.0, 0.0, 0.0, 0.0},
	    []float64{7.0, 6.0, 5.0, 4.0, 6.0, 0.0, 0.0},
	    []float64{7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 0.0},
	    []float64{7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0}}
	   A := matrix.FloatMatrixFromTable(Ldata, matrix.RowOrder)
	   N := A.Rows()
	*/
	N := 7
	nb := 0

	A0 := matrix.FloatUniform(N, N)
	A := matrix.FloatZeros(N, N)
	Mult(A, A0, A0, 1.0, 1.0, TRANSB)

	B := matrix.FloatNormal(A.Rows(), 2)
	B0 := B.Copy()
	B1 := B.Copy()
	Mult(B0, A, B, 1.0, 0.0, NOTRANS)
	_, _, _ = B0, B1, A0

	ipiv := make([]int, N, N)
	L, _ := DecomposeLDL(A.Copy(), nil, ipiv, LOWER, 0)
	//t.Logf("unblk: ipiv = %v\n", ipiv)
	//t.Logf("unblk: L\n%v\n", L)

	ApplyRowPivots(B, ipiv, FORWARD)
	MultTrm(B, L, 1.0, LOWER|UNIT|TRANSA)
	MultDiag(B, L, LEFT)
	MultTrm(B, L, 1.0, LOWER|UNIT)
	ApplyRowPivots(B0, ipiv, FORWARD)
	t.Logf(" unblk: L*D*L.T %d pivots: ||A*B - L*D*L.T*B||_1: %e\n",
		NumPivots(ipiv), NormP(B.Minus(B0), NORM_ONE))
	t.Logf("pivots: %v\n", ipiv)

	nb = 4
	w := matrix.FloatWithValue(A.Rows(), nb, 1.0)
	L, _ = DecomposeLDL(A.Copy(), w, ipiv, LOWER, nb)
	//t.Logf("blk: ipiv = %v\n", ipiv)
	//t.Logf("blk: L\n%v\n", L)

	// B2 = A*B1 == A*B
	B2 := B1.Copy()
	Mult(B2, A, B1, 1.0, 0.0, NOTRANS)

	ApplyRowPivots(B1, ipiv, FORWARD)
	MultTrm(B1, L, 1.0, LOWER|UNIT|TRANSA)
	MultDiag(B1, L, LEFT)
	MultTrm(B1, L, 1.0, LOWER|UNIT)
	ApplyRowPivots(B2, ipiv, FORWARD)
	t.Logf("   blk: L*D*L.T %d pivots: ||A*B - L*D*L.T*B||_1: %e\n",
		NumPivots(ipiv), NormP(B2.Minus(B1), NORM_ONE))
	t.Logf("pivots: %v\n", ipiv)
}

func _TestLDLupper(t *testing.T) {
	N := 8
	nb := 0

	A0 := matrix.FloatUniform(N, N)
	A := matrix.FloatZeros(N, N)
	Mult(A, A0, A0, 1.0, 1.0, TRANSB)

	B := matrix.FloatNormal(A.Rows(), 2)
	B0 := B.Copy()
	B1 := B.Copy()
	Mult(B0, A, B, 1.0, 0.0, NOTRANS)

	ipiv := make([]int, N, N)
	U, _ := DecomposeLDL(A.Copy(), nil, ipiv, UPPER, 0)

	ApplyRowPivots(B, ipiv, BACKWARD)
	MultTrm(B, U, 1.0, UPPER|UNIT|TRANSA)
	MultDiag(B, U, LEFT)
	MultTrm(B, U, 1.0, UPPER|UNIT)
	ApplyRowPivots(B0, ipiv, BACKWARD)
	t.Logf(" unblk: U*D*U.T %d pivots: ||A*B - U*D*U.T*B||_1: %e\n",
		NumPivots(ipiv), NormP(B.Minus(B0), NORM_ONE))
	t.Logf("pivots: %v\n", ipiv)

	nb = 4
	w := matrix.FloatZeros(A.Rows(), nb)
	U, _ = DecomposeLDL(A.Copy(), w, ipiv, UPPER, nb)
	// B2 = A*B1 == A*B
	B2 := B1.Copy()
	Mult(B2, A, B1, 1.0, 0.0, NOTRANS)

	ApplyRowPivots(B1, ipiv, BACKWARD)
	MultTrm(B1, U, 1.0, UPPER|UNIT|TRANSA)
	MultDiag(B1, U, LEFT)
	MultTrm(B1, U, 1.0, UPPER|UNIT)
	ApplyRowPivots(B2, ipiv, BACKWARD)
	t.Logf("   blk: U*D*U.T %d pivots: ||A*B - U*D*U.T*B||_1: %e\n",
		NumPivots(ipiv), NormP(B2.Minus(B1), NORM_ONE))
	t.Logf("pivots: %v\n", ipiv)
}

func TestLDLSolve(t *testing.T) {
	N := 7
	nb := 4
	K := 6
	A0 := matrix.FloatUniform(N, N)
	A := matrix.FloatZeros(N, N)
	Mult(A, A0, A0, 1.0, 1.0, TRANSB)

	X0 := matrix.FloatNormal(A.Rows(), K)
	B0 := matrix.FloatZeros(X0.Size())
	Mult(B0, A, X0, 1.0, 0.0, NOTRANS)
	B := B0.Copy()
	_ = B

	w := matrix.FloatZeros(A.Rows(), nb)
	ipiv := make([]int, N, N)
	L, _ := DecomposeLDL(A.Copy(), w, ipiv, LOWER, nb)

	ipiv0 := make([]int, N, N)
	L0, _ := DecomposeLDL(A.Copy(), w, ipiv0, LOWER, 0)

	t.Logf("L:\n%v\n", L)
	t.Logf("ipiv: %v\n", ipiv)
	t.Logf("L0:\n%v\n", L0)
	t.Logf("ipiv0: %v\n", ipiv0)
	t.Logf("L == L0: %v", L.AllClose(L0))

	// B = A*X0; solve and B should be X0
	/*
	   SolveLDL(B, L, ipiv, LOWER)
	   B.Minus(X0)
	   t.Logf("L*D*L.T: ||A*X - B||_1: %e\n", NormP(B, NORM_ONE))
	   t.Logf("ipiv: %v\n", ipiv)
	*/

	/*
	   B = B0.Copy()
	   U, _ := DecomposeLDL(A.Copy(), w, ipiv, UPPER, nb)
	   // B = A*X0; solve and B should be X0
	   SolveLDL(B, U, ipiv, UPPER)
	   B.Minus(X0)
	   t.Logf("U*D*U.T: ||A*X - B||_1: %e\n", NormP(B, NORM_ONE))
	*/
}

// Local Variables:
// tab-width: 4
// indent-tabs-mode: nil
// End:
