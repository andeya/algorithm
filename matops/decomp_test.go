// Copyright (c) Harri Rautila, 2012,2013

// This file is part of github.com/hrautila/matops package. It is free software,
// distributed under the terms of GNU Lesser General Public License Version 3, or
// any later version. See the COPYING tile included in this archive.

package matops

import (
	"github.com/henrylee2cn/algorithm/linalg"
	"github.com/henrylee2cn/algorithm/linalg/blas"
	"github.com/henrylee2cn/algorithm/linalg/lapack"
	"github.com/henrylee2cn/algorithm/matops/calgo"
	"github.com/henrylee2cn/algorithm/matrix"
	"testing"
	//"math"
)

func _TestLU2x2NoPiv(t *testing.T) {
	Adata2 := [][]float64{
		[]float64{4.0, 3.0},
		[]float64{6.0, 3.0}}

	A := matrix.FloatMatrixFromTable(Adata2, matrix.RowOrder)
	DecomposeBlockSize(0)
	DecomposeLUnoPiv(A, 0)
	t.Logf("A\n%v\n", A)
	Ld := TriLU(A.Copy())
	Ud := TriU(A)
	t.Logf("L*U\n%v\n", matrix.Times(Ld, Ud))
}

func _TestLU3x3NoPiv(t *testing.T) {
	Adata2 := [][]float64{
		[]float64{4.0, 2.0, 2.0},
		[]float64{6.0, 4.0, 2.0},
		[]float64{4.0, 6.0, 1.0},
	}

	A := matrix.FloatMatrixFromTable(Adata2, matrix.RowOrder)
	A0 := A.Copy()
	DecomposeBlockSize(0)
	DecomposeLUnoPiv(A, 0)
	t.Logf("A\n%v\n", A)
	Ld := TriLU(A.Copy())
	Ud := TriU(A.Copy())
	t.Logf("A == L*U: %v\n", A0.AllClose(matrix.Times(Ld, Ud)))
}

func _TestUnblkLUnoPiv(t *testing.T) {
	N := 6
	L := matrix.FloatUniformSymmetric(N, matrix.Lower)
	U := matrix.FloatUniformSymmetric(N, matrix.Upper)
	// Set L diagonal to 1.0
	L.Diag().SetIndexes(1.0)

	A := matrix.Times(L, U)
	t.Logf("A\n%v\n", A)
	DecomposeBlockSize(0)
	R, _ := DecomposeLUnoPiv(A.Copy(), 0)
	Ld := TriLU(R.Copy())
	Ud := TriU(R)
	t.Logf("A == L*U: %v\n", A.AllClose(matrix.Times(Ld, Ud)))
}

func _TestBlkLUnoPiv(t *testing.T) {
	N := 10
	nb := 4
	L := matrix.FloatUniformSymmetric(N, matrix.Lower)
	U := matrix.FloatUniformSymmetric(N, matrix.Upper)
	// Set L diagonal to 1.0
	L.Diag().SetIndexes(1.0)

	A := matrix.Times(L, U)
	t.Logf("A\n%v\n", A)
	DecomposeBlockSize(nb)
	R, _ := DecomposeLUnoPiv(A.Copy(), nb)
	Ld := TriLU(R.Copy())
	Ud := TriU(R)
	t.Logf("A == L*U: %v\n", A.AllClose(matrix.Times(Ld, Ud)))
}

func _TestLU3x3Piv(t *testing.T) {
	Adata2 := [][]float64{
		[]float64{3.0, 2.0, 2.0},
		[]float64{6.0, 4.0, 1.0},
		[]float64{4.0, 6.0, 3.0},
	}
	A := matrix.FloatMatrixFromTable(Adata2, matrix.RowOrder)
	piv := make([]int, A.Rows())
	piv0 := make([]int32, A.Rows())
	A0 := A.Copy()
	t.Logf("start A\n%v\n", A)
	DecomposeBlockSize(0)
	DecomposeLU(A, piv, 0)
	Ld := TriLU(A.Copy())
	Ud := TriU(A.Copy())
	t.Logf("A\n%v\n", A)
	t.Logf("Ld:\n%v\n", Ld)
	t.Logf("Ud:\n%v\n", Ud)
	t.Logf("piv: %v\n", piv)
	t.Logf("result:\n%v\n", matrix.Times(Ld, Ud))
	//t.Logf("A == L*U: %v\n", A0.AllClose(matrix.Times(Ld, Ud)))
	lapack.Getrf(A0, piv0)
	t.Logf("lapack result: piv0 %v\n%v\n", piv0, A0)
	t.Logf("A == A0: %v\n", A0.AllClose(A))
}

func _TestLU3x4Piv(t *testing.T) {
	Adata2 := [][]float64{
		[]float64{3.0, 2.0, 2.0, 1.0},
		[]float64{6.0, 4.0, 1.0, 2.0},
		[]float64{4.0, 6.0, 3.0, 3.0},
	}
	A := matrix.FloatMatrixFromTable(Adata2, matrix.RowOrder)
	piv := make([]int, A.Rows())
	piv0 := make([]int32, A.Rows())
	A0 := A.Copy()
	t.Logf("start A\n%v\n", A)
	DecomposeBlockSize(0)
	DecomposeLU(A, piv, 0)
	t.Logf("piv: %v\n", piv)
	lapack.Getrf(A0, piv0)
	t.Logf("lapack result: piv0 %v\n%v\n", piv0, A0)
	t.Logf("A == A0: %v\n", A0.AllClose(A))
}

func _TestBlkLUPiv(t *testing.T) {
	N := 10
	nb := 4
	L := matrix.FloatUniformSymmetric(N, matrix.Lower)
	U := matrix.FloatUniformSymmetric(N, matrix.Upper)
	// Set L diagonal to 1.0
	L.Diag().SetIndexes(1.0)

	A := matrix.Times(L, U)
	A0 := A.Copy()
	piv := make([]int, N, N)
	DecomposeBlockSize(nb)
	R, _ := DecomposeLU(A.Copy(), piv, 0)
	t.Logf("piv: %v\n", piv)

	piv0 := make([]int32, N, N)
	lapack.Getrf(A0, piv0)
	t.Logf("lapack result: piv0 %v\n", piv0)
	t.Logf("R == A0: %v\n", A0.AllClose(R))
}

func _TestCHOL3x3(t *testing.T) {
	Ldata2 := [][]float64{
		[]float64{3.0, 0.0, 0.0},
		[]float64{6.0, 4.0, 0.0},
		[]float64{4.0, 6.0, 3.0},
	}
	L := matrix.FloatMatrixFromTable(Ldata2, matrix.RowOrder)
	A := matrix.Times(L, L.Transpose())
	DecomposeBlockSize(0)
	DecomposeCHOL(A, LOWER, 0)
	Ld := TriL(A.Copy())
	t.Logf("Ld:\n%v\n", Ld)
	t.Logf("result L == Ld: %v\n", L.AllClose(Ld))
}

func _TestCHOLUpLo(t *testing.T) {
	N := 10
	L := matrix.FloatUniformSymmetric(N, matrix.Lower)
	A := matrix.Times(L, L.Transpose())
	DecomposeCHOL(A, LOWER, 0)
	Ld := TriL(A)
	t.Logf("result L == Ld: %v\n", L.AllClose(Ld))

	U := matrix.FloatUniformSymmetric(N, matrix.Upper)
	A = matrix.Times(U.Transpose(), U)
	DecomposeCHOL(A, UPPER, 0)
	Ud := TriU(A)
	t.Logf("result U == Ud: %v\n", U.AllClose(Ud))
}

func _TestBlkCHOLUpLo(t *testing.T) {
	N := 10
	nb := 4
	L := matrix.FloatUniformSymmetric(N, matrix.Lower)
	A := matrix.Times(L, L.Transpose())
	DecomposeCHOL(A, LOWER, nb)
	Ld := TriL(A)
	ok := L.AllClose(Ld)
	t.Logf("result L == Ld: %v\n", ok)
	if !ok {
		t.Logf("L:\n%v\n", L)
		t.Logf("Ld:\n%v\n", Ld)
	}

	U := matrix.FloatUniformSymmetric(N, matrix.Upper)
	A = matrix.Times(U.Transpose(), U)
	DecomposeCHOL(A, UPPER, nb)
	Ud := TriU(A)
	ok = U.AllClose(Ud)
	t.Logf("result U == Ud: %v\n", ok)

	//lapack.Potrf(A0)
	//t.Logf("lapack result:\n%v\n", A0)
	//t.Logf("A == A0: %v\n", A0.AllClose(A))
}

func TestQRSmal(t *testing.T) {
	data := [][]float64{
		[]float64{12.0, -51.0, 4.0},
		[]float64{6.0, 167.0, -68.0},
		[]float64{-4.0, 24.0, -41.0}}

	A := matrix.FloatMatrixFromTable(data, matrix.RowOrder)
	T := matrix.FloatZeros(A.Cols(), A.Cols())
	T0 := T.Copy()

	M := A.Rows()
	//N := A.Cols()
	Tau := matrix.FloatZeros(M, 1)
	X, _ := DecomposeQR(A.Copy(), Tau, nil, 0)
	t.Logf("A\n%v\n", A)
	t.Logf("X\n%v\n", X)
	t.Logf("Tau\n%v\n", Tau)

	Tau0 := matrix.FloatZeros(M, 1)
	lapack.Geqrf(A, Tau0)
	t.Logf("lapack X\n%v\n", A)
	t.Logf("lapack Tau\n%v\n", Tau0)

	unblkQRBlockReflector(X, Tau, T)
	t.Logf("T:\n%v\n", T)

	V := TriLU(X.Copy())
	lapack.LarftFloat(V, Tau, T0)
	t.Logf("T0:\n%v\n", T0)

}

func TestQR(t *testing.T) {
	M := 6
	N := 5

	A := matrix.FloatUniform(M, N)
	Tau := matrix.FloatZeros(M, 1)

	T := matrix.FloatZeros(A.Cols(), A.Cols())
	T0 := T.Copy()

	X, _ := DecomposeQR(A.Copy(), Tau, nil, 0)

	Tau0 := matrix.FloatZeros(M, 1)
	A0 := A.Copy()
	lapack.Geqrf(A0, Tau0)
	ok := X.AllClose(A0)
	okt := Tau.AllClose(Tau0)
	t.Logf("lapack QR == DecomposeQR: %v\n", ok && okt)
	if !ok || !okt {
		t.Logf("A0: %d, %d, %d\n", A0.Rows(), A0.Cols(), A0.LeadingIndex())
		t.Logf("A\n%v\n", A)
		t.Logf("X\n%v\n", X)
		t.Logf("Tau\n%v\n", Tau)
		t.Logf("lapack X\n%v\n", A0)
		t.Logf("lapack Tau\n%v\n", Tau0)
	}

	// build block reflectors
	unblkQRBlockReflector(X, Tau, T)
	V := TriLU(X.Copy())
	lapack.LarftFloat(V, Tau, T0)

	ok = T0.AllClose(T)
	t.Logf("lapack.dlarft == QRBlockReflector: %v\n", ok)
	if !ok {
		t.Logf("T:\n%v\n", T)
		t.Logf("lapack T0:\n%v\n", T0)
	}
}

func TestQRT(t *testing.T) {
	M := 6
	N := 5

	var Tau matrix.FloatMatrix
	A := matrix.FloatUniform(M, N)
	T := matrix.FloatZeros(A.Cols(), A.Cols())
	T0 := T.Copy()

	X, _ := DecomposeQRT(A.Copy(), T, nil, 0)
	Tau.DiagOf(T)

	Tau0 := matrix.FloatZeros(M, 1)
	A0 := A.Copy()
	lapack.Geqrf(A0, Tau0)
	ok := X.AllClose(A0)
	okt := Tau.AllClose(Tau0)
	t.Logf("lapack QR == DecomposeQR: %v\n", ok && okt)
	if !ok || !okt {
		t.Logf("A0: %d, %d, %d\n", A0.Rows(), A0.Cols(), A0.LeadingIndex())
		t.Logf("A\n%v\n", A)
		t.Logf("X\n%v\n", X)
		t.Logf("Tau\n%v\n", &Tau)
		t.Logf("lapack X\n%v\n", A0)
		t.Logf("lapack Tau\n%v\n", Tau0)
	}

	// build block reflectors
	//unblkQRBlockReflector(X, Tau, T)
	V := TriLU(A0.Copy())
	lapack.LarftFloat(V, Tau0, T0)

	ok = T0.AllClose(T)
	t.Logf("lapack.dlarft == QRBlockReflector: %v\n", ok)
	if !ok {
		t.Logf("T:\n%v\n", T)
		t.Logf("lapack T0:\n%v\n", T0)
	}
}

func TestQRBlk(t *testing.T) {
	M := 8
	N := 6
	nb := 2

	A := matrix.FloatUniform(M, N)
	Tau := matrix.FloatZeros(M, 1)
	Tz := matrix.FloatZeros(N, N)
	Tx := matrix.FloatZeros(N, N)

	DecomposeBlockSize(0)
	Z, _ := DecomposeQRT(A.Copy(), Tz, nil, 0)
	_ = Z

	DecomposeBlockSize(nb)
	//X, _ := DecomposeQR(A.Copy(), Tau, nil)
	X, _ := DecomposeQRT(A.Copy(), Tx, nil, nb)

	Tau0 := matrix.FloatZeros(M, 1)
	A0 := A.Copy()
	lapack.Geqrf(A0, Tau0)
	ok := X.AllClose(A0)
	var dx matrix.FloatMatrix
	dx.DiagOf(Tx)
	//okt := Tau.AllClose(Tau0)
	okt := true
	_ = Tau
	t.Logf("lapack QR == DecomposeQR: %v\n", ok && okt)
	if !ok || !okt || true {
		t.Logf("A0: %d, %d, %d\n", A0.Rows(), A0.Cols(), A0.LeadingIndex())
		t.Logf("A\n%v\n", A)
		t.Logf("X\n%v\n", X)
		t.Logf("Tz\n%v\n", Tz)
		t.Logf("Tx\n%v\n", Tx)
		t.Logf("Tau\n%v\n", &dx)
		t.Logf("lapack X\n%v\n", A0)
		t.Logf("lapack Tau\n%v\n", Tau0)
	}
}

func _TestBlkQRUT(t *testing.T) {
	data := [][]float64{
		[]float64{12.0, -51.0, 4.0},
		[]float64{6.0, 167.0, -68.0},
		[]float64{-4.0, 24.0, -41.0}}

	A := matrix.FloatMatrixFromTable(data, matrix.RowOrder)

	M := A.Rows()
	N := A.Cols()
	var d matrix.FloatMatrix
	//nb := 0
	T := matrix.FloatZeros(N, N)
	Q := matrix.FloatZeros(M, M)
	X, _ := DecomposeQRUT(A.Copy(), T)
	t.Logf("A\n%v\n", A)
	Y := TriLU(X.Copy())
	R := TriU(X.Copy())
	t.Logf("R\n%v\n", R)
	t.Logf("Y\n%v\n", Y)
	t.Logf("T\n%v\n", T)
	_ = Q
	_ = d
}

func update(t *testing.T, Y1, Y2, C1, C2, T, W *matrix.FloatMatrix) {
	if W.Rows() != C1.Cols() {
		panic("W.Rows != C1.Cols")
	}
	// W = C1.T
	ScalePlus(W, C1, 0.0, 1.0, TRANSB)
	//fmt.Printf("W = C1.T:\n%v\n", W)
	// W = C1.T*Y1
	//MultTrm(W, Y1, 1.0, LOWER|UNIT|RIGHT)
	Wr := W.FloatArray()
	Y1r := Y1.FloatArray()
	ldW := W.LeadingIndex()
	ldY := Y1.LeadingIndex()
	calgo.DTrmmUnblk(Wr, Y1r, 1.0, calgo.Flags(LOWER|UNIT|RIGHT),
		ldW, ldY, Y1.Cols(), 0, W.Rows(), 0)
	t.Logf("W = C1.T*Y1:\n%v\n", W)
	// W = W + C2.T*Y2
	Mult(W, C2, Y2, 1.0, 1.0, TRANSA)
	t.Logf("W = W + C2.T*Y2:\n%v\n", W)

	// --- here: W == C.T*Y ---
	// W = W*T
	MultTrm(W, T, 1.0, UPPER|RIGHT)
	t.Logf("W = C.T*Y*T:\n%v\n", W)

	// --- here: W == C.T*Y*T ---
	// C2 = C2 - Y2*W.T
	Mult(C2, Y2, W, -1.0, 1.0, TRANSB)
	t.Logf("C2 = C2 - Y2*W.T:\n%v\n", C2)
	//  W = Y1*W.T ==> W.T = W*Y1.T
	MultTrm(W, Y1, 1.0, LOWER|UNIT|TRANSA|RIGHT)
	t.Logf("W.T = W*Y1.T:\n%v\n", W)

	// C1 = C1 - W.T
	ScalePlus(C1, W, 1.0, -1.0, TRANSB)
	//fmt.Printf("C1 = C1 - W.T:\n%v\n", C1)

	// --- here: C = (I - Y*T*Y.T).T * C ---
}

func updateBlas(t *testing.T, Y1, Y2, C1, C2, T, W *matrix.FloatMatrix) {
	if W.Rows() != C1.Cols() {
		panic("W.Rows != C1.Cols")
	}
	// W = C1.T
	ScalePlus(W, C1, 0.0, 1.0, TRANSB)
	//fmt.Printf("W = C1.T:\n%v\n", W)
	// W = C1.T*Y1
	blas.TrmmFloat(Y1, W, 1.0, linalg.OptLower, linalg.OptUnit, linalg.OptRight)
	t.Logf("W = C1.T*Y1:\n%v\n", W)
	// W = W + C2.T*Y2
	blas.GemmFloat(C2, Y2, W, 1.0, 1.0, linalg.OptTransA)
	t.Logf("W = W + C2.T*Y2:\n%v\n", W)

	// --- here: W == C.T*Y ---
	// W = W*T
	blas.TrmmFloat(T, W, 1.0, linalg.OptUpper, linalg.OptRight)
	t.Logf("W = C.T*Y*T:\n%v\n", W)

	// --- here: W == C.T*Y*T ---
	// C2 = C2 - Y2*W.T
	blas.GemmFloat(Y2, W, C2, -1, 1.0, linalg.OptTransB)
	t.Logf("C2 = C2 - Y2*W.T:\n%v\n", C2)
	//  W = Y1*W.T ==> W.T = W*Y1.T
	blas.TrmmFloat(Y1, W, 1.0, linalg.OptLower,
		linalg.OptUnit, linalg.OptRight, linalg.OptTrans)
	t.Logf("W.T = W*Y1.T:\n%v\n", W)

	// C1 = C1 - W.T
	ScalePlus(C1, W, 1.0, -1.0, TRANSB)
	//fmt.Printf("C1 = C1 - W.T:\n%v\n", C1)

	// --- here: C = (I - Y*T*Y.T).T * C ---
}

func _TestUpdate(t *testing.T) {
	Tdata := [][]float64{
		[]float64{1.37e+00, 3.77e-01},
		[]float64{0.00e+00, 1.39e+00}}

	C1data := [][]float64{
		[]float64{2.54e-01, 9.77e-01, 8.01e-01, 9.08e-02},
		[]float64{2.82e-01, 7.43e-02, 7.30e-01, 4.93e-01}}

	C2data := [][]float64{
		[]float64{7.89e-01, 2.22e-01, 1.83e-01, 9.27e-01},
		[]float64{3.62e-01, 6.81e-01, 4.28e-01, 9.55e-01},
		[]float64{8.81e-01, 2.42e-01, 8.97e-01, 3.48e-01},
		[]float64{2.97e-01, 3.12e-01, 6.83e-01, 6.91e-01},
		[]float64{8.94e-01, 9.33e-01, 9.79e-01, 7.11e-01},
		[]float64{9.75e-02, 7.42e-01, 9.22e-01, 5.64e-01}}

	Y1data := [][]float64{
		[]float64{-1.41e+00, -1.11e+00},
		[]float64{1.46e-02, -7.04e-01}}

	Y2data := [][]float64{
		[]float64{8.19e-02, 2.62e-02},
		[]float64{3.14e-01, -2.57e-02},
		[]float64{5.05e-01, -3.73e-01},
		[]float64{4.11e-02, 2.09e-01},
		[]float64{3.08e-01, -1.34e-01},
		[]float64{3.06e-02, 4.86e-01}}

	Wdata := [][]float64{
		[]float64{1.62e+00, 3.81e-01},
		[]float64{2.28e+00, 1.01e+00},
		[]float64{2.42e+00, 1.84e+00},
		[]float64{1.24e+00, 1.30e+00}}

	T := matrix.FloatMatrixFromTable(Tdata, matrix.RowOrder)
	C1 := matrix.FloatMatrixFromTable(C1data, matrix.RowOrder)
	Y1 := matrix.FloatMatrixFromTable(Y1data, matrix.RowOrder)
	C2 := matrix.FloatMatrixFromTable(C2data, matrix.RowOrder)
	Y2 := matrix.FloatMatrixFromTable(Y2data, matrix.RowOrder)
	W := matrix.FloatMatrixFromTable(Wdata, matrix.RowOrder)
	_, _, _ = C2, Y2, W
	_, _, _ = T, C1, Y1

	C1b := C1.Copy()
	C2b := C2.Copy()
	update(t, Y1, Y2, C1, C2, T, W)
	t.Logf("C1:\n%v\n", C1)
	updateBlas(t, Y1, Y2, C1b, C2b, T, W)
	t.Logf("C1b:\n%v\n", C1b)

}

// Local Variables:
// tab-width: 4
// indent-tabs-mode: nil
// End:
