// Copyright (c) Harri Rautila, 2012,2013

// This file is part of github.com/hrautila/matops package. It is free software,
// distributed under the terms of GNU Lesser General Public License Version 3, or
// any later version. See the COPYING tile included in this archive.

package calgo

import (
	"github.com/henrylee2cn/algorithm/linalg"
	"github.com/henrylee2cn/algorithm/linalg/blas"
	"github.com/henrylee2cn/algorithm/matrix"
	"testing"
	//"math/rand"
	//"math"
	//"time"
)

func _TestMultSmall(t *testing.T) {
	bM := 6
	bN := 6
	bP := 6
	D := matrix.FloatNormal(bM, bP)
	E := matrix.FloatNormal(bP, bN)
	C0 := matrix.FloatWithValue(bM, bN, 1.0)
	C1 := C0.Copy()

	Dr := D.FloatArray()
	Er := E.FloatArray()
	C1r := C1.FloatArray()

	blas.GemmFloat(D, E, C0, 1.0, 1.0)
	t.Logf("blas: C=D*E\n%v\n", C0)

	DMult(C1r, Dr, Er, 1.0, 1.0, NOTRANS, bM, bM, bP, bP, 0, bN, 0, bM, 4, 4, 4)
	t.Logf("C0 == C1: %v\n", C0.AllClose(C1))
	t.Logf("C1: C1=D*E\n%v\n", C1)
}

func _TestMultBig(t *testing.T) {
	bM := 100*M + 3
	bN := 100*N + 3
	bP := 100*P + 3
	D := matrix.FloatNormal(bM, bP)
	E := matrix.FloatNormal(bP, bN)
	C0 := matrix.FloatZeros(bM, bN)
	C1 := matrix.FloatZeros(bM, bN)

	Dr := D.FloatArray()
	Er := E.FloatArray()
	C1r := C1.FloatArray()

	blas.GemmFloat(D, E, C0, 1.0, 1.0)
	//t.Logf("blas: C=D*E\n%v\n", C0)

	DMult(C1r, Dr, Er, 1.0, 1.0, NOTRANS, bM, bM, bP, bP, 0, bN, 0, bM, 32, 32, 32)
	res := C0.AllClose(C1)
	t.Logf("C0 == C1: %v\n", res)
}

func _TestMultTransASmall(t *testing.T) {
	bM := 7
	bN := 7
	bP := 7
	D := matrix.FloatNormal(bM, bP)
	E := matrix.FloatNormal(bP, bN)
	C0 := matrix.FloatWithValue(bM, bN, 0.0)
	C1 := C0.Copy()
	Dt := D.Transpose()

	Dr := Dt.FloatArray()
	Er := E.FloatArray()
	C1r := C1.FloatArray()
	blas.GemmFloat(Dt, E, C0, 1.0, 1.0, linalg.OptTransA)
	t.Logf("blas: C=D*E\n%v\n", C0)

	DMult(C1r, Dr, Er, 1.0, 1.0, TRANSA, bM, bM, bP, bP, 0, bN, 0, bM, 4, 4, 4)
	t.Logf("C0 == C1: %v\n", C0.AllClose(C1))
	t.Logf("C1: C1=D*E\n%v\n", C1)
}

func _TestMultTransABig(t *testing.T) {
	bM := 100*M + 3
	bN := 100*N + 3
	bP := 100*P + 3
	D := matrix.FloatNormal(bM, bP)
	E := matrix.FloatNormal(bP, bN)
	C0 := matrix.FloatZeros(bM, bN)
	C1 := matrix.FloatZeros(bM, bN)
	Dt := D.Transpose()

	Dr := Dt.FloatArray()
	Er := E.FloatArray()
	C1r := C1.FloatArray()

	blas.GemmFloat(Dt, E, C0, 1.0, 1.0, linalg.OptTransA)

	DMult(C1r, Dr, Er, 1.0, 1.0, TRANSA, bM, bM, bP, bP, 0, bN, 0, bM, 32, 32, 32)
	t.Logf("C0 == C1: %v\n", C0.AllClose(C1))
}

func _TestMultTransBSmall(t *testing.T) {
	bM := 7
	bN := 7
	bP := 7
	D := matrix.FloatNormal(bM, bP)
	E := matrix.FloatNormal(bP, bN)
	C0 := matrix.FloatWithValue(bP, bN, 1.0)
	C1 := C0.Copy()
	Et := E.Transpose()

	Dr := D.FloatArray()
	Er := Et.FloatArray()
	C1r := C1.FloatArray()

	blas.GemmFloat(D, Et, C0, 1.0, 1.0, linalg.OptTransB)
	t.Logf("blas: C=D*E.T\n%v\n", C0)

	DMult(C1r, Dr, Er, 1.0, 1.0, TRANSB, bM, bM, bP, bP, 0, bN, 0, bM, 4, 4, 4)
	t.Logf("C0 == C1: %v\n", C0.AllClose(C1))
	t.Logf("C1: C1=D*E.T\n%v\n", C1)
}

func _TestMultTransBBig(t *testing.T) {
	bM := 100*M + 3
	bN := 100*N + 3
	bP := 100*P + 3
	D := matrix.FloatNormal(bM, bP)
	E := matrix.FloatNormal(bP, bN)
	C0 := matrix.FloatZeros(bM, bN)
	C1 := matrix.FloatZeros(bM, bN)
	Et := E.Transpose()

	Dr := D.FloatArray()
	Er := Et.FloatArray()
	C1r := C1.FloatArray()

	blas.GemmFloat(D, Et, C0, 1.0, 1.0, linalg.OptTransB)

	DMult(C1r, Dr, Er, 1.0, 1.0, TRANSB, bM, bM, bP, bP, 0, bN, 0, bM, 32, 32, 32)
	t.Logf("C0 == C1: %v\n", C0.AllClose(C1))
}

func _TestMultTransABSmall(t *testing.T) {
	bM := 7
	bN := 7
	bP := 7
	D := matrix.FloatNormal(bM, bP)
	E := matrix.FloatNormal(bP, bN)
	C0 := matrix.FloatZeros(bM, bN)
	C1 := matrix.FloatZeros(bM, bN)
	Dt := D.Transpose()
	Et := E.Transpose()

	Dr := Dt.FloatArray()
	Er := Et.FloatArray()
	C1r := C1.FloatArray()

	blas.GemmFloat(Dt, Et, C0, 1.0, 1.0, linalg.OptTransA, linalg.OptTransB)
	t.Logf("blas: C=D.T*E.T\n%v\n", C0)

	DMult(C1r, Dr, Er, 1.0, 1.0, TRANSA|TRANSB, bM, bM, bP, bP, 0, bN, 0, bM, 4, 4, 4)
	t.Logf("C0 == C1: %v\n", C0.AllClose(C1))
	t.Logf("C1: C1=D.T*E.T\n%v\n", C1)
}

func _TestMultTransABBig(t *testing.T) {
	bM := 100 * M
	bN := 100 * N
	bP := 100 * P
	D := matrix.FloatNormal(bM, bP)
	E := matrix.FloatNormal(bP, bN)
	C0 := matrix.FloatZeros(bM, bN)
	C1 := matrix.FloatZeros(bM, bN)
	Dt := D.Transpose()
	Et := E.Transpose()

	Dr := Dt.FloatArray()
	Er := Et.FloatArray()
	C1r := C1.FloatArray()

	blas.GemmFloat(Dt, Et, C0, 1.0, 1.0, linalg.OptTransA, linalg.OptTransB)

	DMult(C1r, Dr, Er, 1.0, 1.0, TRANSA|TRANSB, bM, bM, bP, bP, 0, bN, 0, bM, 32, 32, 32)
	t.Logf("C0 == C1: %v\n", C0.AllClose(C1))
}

func _TestMultMVSmall(t *testing.T) {
	bM := 5
	bN := 4
	A := matrix.FloatNormal(bM, bN)
	X := matrix.FloatVector([]float64{1.0, 2.0, 3.0, 4.0})
	Y1 := matrix.FloatZeros(bM, 1)
	Y0 := matrix.FloatZeros(bM, 1)

	Ar := A.FloatArray()
	Xr := X.FloatArray()
	Y1r := Y1.FloatArray()

	blas.GemvFloat(A, X, Y0, 1.0, 1.0)

	DMultMV(Y1r, Ar, Xr, 1.0, 1.0, NOTRANS, 1, A.LeadingIndex(), 1, 0, bN, 0, bM, 4, 4)
	ok := Y0.AllClose(Y1)
	t.Logf("Y0 == Y1: %v\n", ok)
	if !ok {
		t.Logf("blas: Y=A*X\n%v\n", Y0)
		t.Logf("Y1: Y1 = A*X\n%v\n", Y1)
	}
}

func _TestMultMV(t *testing.T) {
	bM := 100 * M
	bN := 100 * N
	A := matrix.FloatNormal(bM, bN)
	X := matrix.FloatNormal(bN, 1)
	Y1 := matrix.FloatZeros(bM, 1)
	Y0 := matrix.FloatZeros(bM, 1)

	Ar := A.FloatArray()
	Xr := X.FloatArray()
	Y1r := Y1.FloatArray()

	blas.GemvFloat(A, X, Y0, 1.0, 1.0)

	DMultMV(Y1r, Ar, Xr, 1.0, 1.0, NOTRANS, 1, A.LeadingIndex(), 1, 0, bN, 0, bM, 32, 32)
	t.Logf("Y0 == Y1: %v\n", Y0.AllClose(Y1))
	/*
	   if ! Y0.AllClose(Y1) {
	       y0 := Y0.SubMatrix(0, 0, 2, 1)
	       y1 := Y1.SubMatrix(0, 0, 2, 1)
	       t.Logf("y0=\n%v\n", y0)
	       t.Logf("y1=\n%v\n", y1)
	   }
	*/
}

func TestMultMVTransASmall(t *testing.T) {
	data6 := [][]float64{
		[]float64{-1.59e+00, 6.56e-02, 2.14e-01, 6.79e-01, 2.93e-01, 5.24e-01},
		[]float64{4.28e-01, 1.57e-01, 3.81e-01, 2.19e-01, 2.97e-01, 2.83e-02},
		[]float64{3.02e-01, 9.70e-02, 3.18e-01, 2.03e-01, 7.53e-01, 1.58e-01},
		[]float64{1.99e-01, 3.01e-01, 4.69e-01, 3.61e-01, 2.07e-01, 6.07e-01},
		[]float64{1.93e-01, 5.15e-01, 2.83e-01, 5.71e-01, 8.65e-01, 9.75e-01},
		[]float64{3.13e-01, 8.14e-01, 2.93e-01, 8.62e-01, 6.97e-01, 7.95e-02}}
	data5 := [][]float64{
		[]float64{1.57e-01, 3.81e-01, 2.19e-01, 2.97e-01, 2.83e-02},
		[]float64{9.70e-02, 3.18e-01, 2.03e-01, 7.53e-01, 1.58e-01},
		[]float64{3.01e-01, 4.69e-01, 3.61e-01, 2.07e-01, 6.07e-01},
		[]float64{5.15e-01, 2.83e-01, 5.71e-01, 8.65e-01, 9.75e-01},
		[]float64{8.14e-01, 2.93e-01, 8.62e-01, 6.97e-01, 7.95e-02}}
	data2 := []float64{4.28e-01, 3.02e-01, 1.99e-01, 1.93e-01, 3.13e-01}

	bM := 5
	bN := 4
	nb := 2
	//A := matrix.FloatNormal(bN, bM)
	//X := matrix.FloatWithValue(bN, 1, 1.0)

	A := matrix.FloatMatrixFromTable(data5, matrix.RowOrder)
	X := matrix.FloatNew(5, 1, data2)
	bM = A.Rows()
	bN = A.Cols()
	Ym := matrix.FloatZeros(3, bM)
	Y1 := matrix.FloatZeros(bM, 1)
	Y0 := matrix.FloatZeros(bM, 1)

	Ar := A.FloatArray()
	Xr := X.FloatArray()
	Y1r := Y1.FloatArray()

	blas.GemvFloat(A, X, Y0, 1.0, 1.0, linalg.OptTrans)

	DMultMV(Y1r, Ar, Xr, 1.0, 1.0, TRANSA, 1, A.LeadingIndex(), 1, 0, bN, 0, bM, nb, nb)
	ok := Y0.AllClose(Y1)
	t.Logf("Y0 == Y1: %v\n", ok)
	if ok || !ok {
		t.Logf("blas: Y=A.T*X\n%v\n", Y0)
		t.Logf("Y1: Y1 = A*X\n%v\n", Y1)
	}

	// zero Y0, Y1
	Y0.Scale(0.0)
	Y1.Scale(0.0)

	// test with matrix view; A is view
	var A0 matrix.FloatMatrix
	A6 := matrix.FloatMatrixFromTable(data6, matrix.RowOrder)
	A0.SubMatrixOf(A6, 1, 1)

	blas.GemvFloat(&A0, X, Y0, 1.0, 1.0, linalg.OptTrans)

	Ar = A0.FloatArray()
	DMultMV(Y1r, Ar, Xr, 1.0, 1.0, TRANSA, 1, A0.LeadingIndex(), 1, 0, bN, 0, bM, nb, nb)
	ok = Y0.AllClose(Y1)
	t.Logf("lda>rows: Y0 == Y1: %v\n", ok)
	if ok || !ok {
		t.Logf("blas: Y=A.T*X\n%v\n", Y0)
		t.Logf("Y1: Y1 = A*X\n%v\n", Y1)
	}

	// Y is view too.
	Y1.SubMatrixOf(Ym, 0, 0, 1, bM)
	Y1r = Y1.FloatArray()
	DMultMV(Y1r, Ar, Xr, 1.0, 1.0, TRANSA, Y1.LeadingIndex(), A0.LeadingIndex(), 1, 0, bN, 0, bM, nb, nb)
	ok = Y0.AllClose(Y1.Transpose())
	t.Logf("Y0 == Y1 row: %v\n", ok)
	t.Logf("row Y1: %v\n", Y1)
}

func _TestMultMVTransA(t *testing.T) {
	bM := 1000 * M
	bN := 1000 * N
	A := matrix.FloatNormal(bN, bM)
	X := matrix.FloatWithValue(bN, 1, 1.0)
	Y1 := matrix.FloatZeros(bM, 1)
	Y0 := matrix.FloatZeros(bM, 1)

	Ar := A.FloatArray()
	Xr := X.FloatArray()
	Y1r := Y1.FloatArray()

	blas.GemvFloat(A, X, Y0, 1.0, 1.0, linalg.OptTrans)

	DMultMV(Y1r, Ar, Xr, 1.0, 1.0, TRANSA, 1, A.LeadingIndex(), 1, 0, bN, 0, bM, 4, 4)
	ok := Y0.AllClose(Y1)
	t.Logf("Y0 == Y1: %v\n", ok)
	if !ok {
		var y1, y0 matrix.FloatMatrix
		Y1.SubMatrix(&y1, 0, 0, 5, 1)
		t.Logf("Y1[0:5]:\n%v\n", y1)
		Y0.SubMatrix(&y0, 0, 0, 5, 1)
		t.Logf("Y0[0:5]:\n%v\n", y0)
	}
}

func _TestMultSymmSmall(t *testing.T) {
	//bM := 5
	bN := 7
	bP := 7
	Adata := [][]float64{
		[]float64{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
		[]float64{0.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0},
		[]float64{0.0, 0.0, 3.0, 3.0, 3.0, 3.0, 3.0},
		[]float64{0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 4.0},
		[]float64{0.0, 0.0, 0.0, 0.0, 5.0, 5.0, 5.0},
		[]float64{0.0, 0.0, 0.0, 0.0, 0.0, 6.0, 6.0},
		[]float64{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.0}}

	A := matrix.FloatMatrixFromTable(Adata, matrix.RowOrder)
	B := matrix.FloatWithValue(bN, bP, 2.0)
	C0 := matrix.FloatZeros(bN, bP)
	C1 := matrix.FloatZeros(bN, bP)

	Ar := A.FloatArray()
	Br := B.FloatArray()
	C1r := C1.FloatArray()

	blas.SymmFloat(A, B, C0, 1.0, 1.0, linalg.OptUpper, linalg.OptRight)

	DMultSymm(C1r, Ar, Br, 1.0, 1.0, UPPER|RIGHT, bN, A.LeadingIndex(), bN, bN, 0, bP, 0, bN, 2, 2, 2)
	ok := C0.AllClose(C1)
	t.Logf("C0 == C1: %v\n", ok)
	if !ok {
		t.Logf("A=\n%v\n", A)
		t.Logf("blas: C=A*B\n%v\n", C0)
		t.Logf("C1: C1 = A*X\n%v\n", C1)
	}
}

func _TestMultSymmLowerSmall(t *testing.T) {
	//bM := 5
	bN := 7
	bP := 7
	Adata := [][]float64{
		[]float64{1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		[]float64{1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		[]float64{1.0, 2.0, 3.0, 0.0, 0.0, 0.0, 0.0},
		[]float64{1.0, 2.0, 3.0, 4.0, 0.0, 0.0, 0.0},
		[]float64{1.0, 2.0, 3.0, 4.0, 5.0, 0.0, 0.0},
		[]float64{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 0.0},
		[]float64{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0}}

	A := matrix.FloatMatrixFromTable(Adata, matrix.RowOrder)
	B := matrix.FloatNormal(bN, bP)
	C0 := matrix.FloatZeros(bN, bP)
	C1 := matrix.FloatZeros(bN, bP)

	Ar := A.FloatArray()
	Br := B.FloatArray()
	C1r := C1.FloatArray()

	blas.SymmFloat(A, B, C0, 1.0, 1.0, linalg.OptLower, linalg.OptRight)

	DMultSymm(C1r, Ar, Br, 1.0, 1.0, LOWER|RIGHT, bN, A.LeadingIndex(), bN,
		bN, 0, bP, 0, bN, 2, 2, 2)
	ok := C0.AllClose(C1)
	t.Logf("C0 == C1: %v\n", ok)
	if !ok {
		t.Logf("A=\n%v\n", A)
		t.Logf("blas: C=A*B\n%v\n", C0)
		t.Logf("C1: C1 = A*X\n%v\n", C1)
	}
}

func _TestMultSymmUpper(t *testing.T) {
	//bM := 5
	bN := 100*N + 3
	bP := 100*P + 3
	A := matrix.FloatNormalSymmetric(bN, matrix.Upper)
	B := matrix.FloatNormal(bN, bP)
	C0 := matrix.FloatZeros(bN, bP)
	C1 := matrix.FloatZeros(bN, bP)

	Ar := A.FloatArray()
	Br := B.FloatArray()
	C1r := C1.FloatArray()

	blas.SymmFloat(A, B, C0, 1.0, 1.0, linalg.OptUpper)

	DMultSymm(C1r, Ar, Br, 1.0, 1.0, UPPER|LEFT, bN, A.LeadingIndex(), bN,
		bN, 0, bP, 0, bN, 32, 32, 32)
	t.Logf("C0 == C1: %v\n", C0.AllClose(C1))
}

func _TestMultSymmLower(t *testing.T) {
	//bM := 5
	bN := 100 * N
	bP := 100 * P
	A := matrix.FloatNormalSymmetric(bN, matrix.Lower)
	B := matrix.FloatNormal(bN, bP)
	C0 := matrix.FloatZeros(bN, bP)
	C1 := matrix.FloatZeros(bN, bP)

	Ar := A.FloatArray()
	Br := B.FloatArray()
	C1r := C1.FloatArray()

	blas.SymmFloat(A, B, C0, 1.0, 1.0, linalg.OptLower)

	DMultSymm(C1r, Ar, Br, 1.0, 1.0, LOWER|LEFT, bN, A.LeadingIndex(), bN,
		bN, 0, bP, 0, bN, 32, 32, 32)
	t.Logf("C0 == C1: %v\n", C0.AllClose(C1))
}

func TestUpdateTrmLower(t *testing.T) {
	//bM := 5
	bN := 8
	bP := 4
	nb := 4
	A := matrix.FloatNormal(bN, bP)
	//B := matrix.FloatNormal(bP, bN)
	B := A.Transpose()
	C0 := matrix.FloatZeros(bN, bN)
	C2 := matrix.FloatZeros(bN, bN)
	C1 := matrix.FloatZeros(bN, bN)

	Ar := A.FloatArray()
	Br := B.FloatArray()
	C1r := C1.FloatArray()
	C0r := C0.FloatArray()
	C2r := C2.FloatArray()

	// no transpose
	DMult(C1r, Ar, Br, 2.0, 1.0, NOTRANS, C1.LeadingIndex(), A.LeadingIndex(), B.LeadingIndex(),
		bP, 0, bN, 0, bN, nb, nb, nb)
	DTrmUpdBlk(C0r, Ar, Br, 2.0, 1.0, LOWER, C0.LeadingIndex(), A.LeadingIndex(), B.LeadingIndex(),
		bP, 0, bN, nb, nb)
	DTrmUpdBlk(C2r, Ar, Br, 2.0, 1.0, UPPER, C2.LeadingIndex(), A.LeadingIndex(), B.LeadingIndex(),
		bP, 0, bN, nb, nb)

	//t.Logf("C1:\n%v\nC0:\n%v\nC2:\n%v\n", C1, C0, C2)
	// C0 == C2.T
	t.Logf("C0 == C2.T: %v\n", C0.AllClose(C2.Transpose()))
	// C1 == C1 - C2 + C0.T
	Cn := matrix.Minus(C1, C2)
	Cn.Plus(C0.Transpose())
	t.Logf("C1 == C1 - C2 + C0.T: %v\n", Cn.AllClose(C1))

	// B == A.T
	DMult(C1r, Ar, Ar, 2.0, 0.0, TRANSB, C1.LeadingIndex(), A.LeadingIndex(), A.LeadingIndex(),
		bP, 0, bN, 0, bN, nb, nb, nb)
	DTrmUpdBlk(C0r, Ar, Ar, 2.0, 0.0, LOWER|TRANSB, C0.LeadingIndex(), A.LeadingIndex(), A.LeadingIndex(),
		bP, 0, bN, nb, nb)
	DTrmUpdBlk(C2r, Ar, Ar, 2.0, 0.0, UPPER|TRANSB, C2.LeadingIndex(), A.LeadingIndex(), A.LeadingIndex(),
		bP, 0, bN, nb, nb)

	//t.Logf("TRANSB:\nC1:\n%v\nC0:\n%v\nC2:\n%v\n", C1, C0, C2)
	// C0 == C2.T
	t.Logf("B.T: C0 == C2.T: %v\n", C0.AllClose(C2.Transpose()))
	// C1 == C1 - C2 + C0.T
	Cn = matrix.Minus(C1, C2)
	Cn.Plus(C0.Transpose())
	t.Logf("B.T: C1 == C1 - C2 + C0.T: %v\n", Cn.AllClose(C1))

	// A == B.T
	DMult(C1r, Br, Br, 2.0, 0.0, TRANSA, C1.LeadingIndex(), B.LeadingIndex(), B.LeadingIndex(),
		bP, 0, bN, 0, bN, nb, nb, nb)
	DTrmUpdBlk(C0r, Br, Br, 2.0, 0.0, LOWER|TRANSA, C0.LeadingIndex(), B.LeadingIndex(), B.LeadingIndex(),
		bP, 0, bN, nb, nb)
	DTrmUpdBlk(C2r, Br, Br, 2.0, 0.0, UPPER|TRANSA, C2.LeadingIndex(), B.LeadingIndex(), B.LeadingIndex(),
		bP, 0, bN, nb, nb)

	//t.Logf("TRANSA:\nC1:\n%v\nC0:\n%v\nC2:\n%v\n", C1, C0, C2)
	// C0 == C2.T
	t.Logf("A.T: C0 == C2.T: %v\n", C0.AllClose(C2.Transpose()))
	// C1 == C1 - C2 + C0.T
	Cn = matrix.Minus(C1, C2)
	Cn.Plus(C0.Transpose())
	t.Logf("A.T: C1 == C1 - C2 + C0.T: %v\n", Cn.AllClose(C1))

	// A == B.T, B == A.T
	DMult(C1r, Br, Ar, 2.0, 0.0, TRANSA|TRANSB, C1.LeadingIndex(), B.LeadingIndex(), A.LeadingIndex(),
		bP, 0, bN, 0, bN, nb, nb, nb)
	DTrmUpdBlk(C0r, Br, Ar, 2.0, 0.0, LOWER|TRANSA|TRANSB, C0.LeadingIndex(), B.LeadingIndex(), A.LeadingIndex(),
		bP, 0, bN, nb, nb)
	DTrmUpdBlk(C2r, Br, Ar, 2.0, 0.0, UPPER|TRANSA|TRANSB, C2.LeadingIndex(), B.LeadingIndex(), A.LeadingIndex(),
		bP, 0, bN, nb, nb)

	//t.Logf("TRANSA|TRANSB:\nC1:\n%v\nC0:\n%v\nC2:\n%v\n", C1, C0, C2)
	// C0 == C2.T
	t.Logf("A.T, B.T: C0 == C2.T: %v\n", C0.AllClose(C2.Transpose()))
	// C1 == C1 - C2 + C0.T
	Cn = matrix.Minus(C1, C2)
	Cn.Plus(C0.Transpose())
	t.Logf("A.T, B.T: C1 == C1 - C2 + C0.T: %v\n", Cn.AllClose(C1))

}

func TestUpdateTrmMV(t *testing.T) {
	//bM := 5
	bN := 8
	//bP := 4
	nb := 4
	X := matrix.FloatNormal(bN, 1)
	//B := matrix.FloatNormal(bP, bN)
	Y := X.Copy()
	C0 := matrix.FloatZeros(bN, bN)
	C2 := matrix.FloatZeros(bN, bN)
	C1 := matrix.FloatZeros(bN, bN)

	Xr := X.FloatArray()
	Yr := Y.FloatArray()
	C1r := C1.FloatArray()
	C0r := C0.FloatArray()
	C2r := C2.FloatArray()

	// no transpose
	DRankMV(C1r, Xr, Yr, 1.0, C1.LeadingIndex(), 1, 1,
		0, bN, 0, bN, nb, nb)
	DTrmUpdMV(C0r, Xr, Yr, 1.0, LOWER, C0.LeadingIndex(), 1, 1,
		0, bN, nb)
	DTrmUpdMV(C2r, Xr, Yr, 1.0, UPPER, C2.LeadingIndex(), 1, 1,
		0, bN, nb)

	t.Logf("C1:\n%v\nC0:\n%v\nC2:\n%v\n", C1, C0, C2)
	// C0 == C2.T
	t.Logf("C0 == C2.T: %v\n", C0.AllClose(C2.Transpose()))
	// C1 == C1 - C2 + C0.T
	Cn := matrix.Minus(C1, C2)
	Cn.Plus(C0.Transpose())
	t.Logf("C1 == C1 - C2 + C0.T: %v\n", Cn.AllClose(C1))

}

// Local Variables:
// tab-width: 4
// indent-tabs-mode: nil
// End:
