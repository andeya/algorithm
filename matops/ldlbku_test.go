// Copyright (c) Harri Rautila, 2012,2013

// This file is part of github.com/hrautila/matops package. It is free software,
// distributed under the terms of GNU Lesser General Public License Version 3, or
// any later version. See the COPYING tile included in this archive.

package matops

import (
	"github.com/henrylee2cn/algorithm/linalg"
	"github.com/henrylee2cn/algorithm/linalg/lapack"
	"github.com/henrylee2cn/algorithm/matrix"
	"testing"
)

func ldlbkDecomposeUpperTest(A *matrix.FloatMatrix, t *testing.T) {

	TriU(A)
	N := A.Cols()
	nb := 0
	W := matrix.FloatZeros(A.Rows(), 5)
	ipiv := make([]int, N, N)
	L, _ := DecomposeBK(A.Copy(), W, ipiv, UPPER, nb)
	t.Logf("L:\n%v\n", L)
	t.Logf("unblked ipiv: %v\n", ipiv)

	ipiv0 := make([]int, N, N)
	nb = 4
	L0, _ := DecomposeBK(A.Copy(), W, ipiv0, UPPER, nb)
	t.Logf("L:\n%v\n", L0)
	t.Logf("blked ipiv: %v\n", ipiv0)

	ipiv2 := make([]int32, N, N)
	lapack.SytrfFloat(A, ipiv2, linalg.OptUpper)
	t.Logf("lapack A:\n%v\n", A)
	t.Logf("lapack ipiv: %v\n", ipiv2)

}

func _TestBKSavedUpper(t *testing.T) {
	Aorig := fromFile("test/saved.dat", t)
	ldlbkDecomposeUpperTest(Aorig, t)
}

func _TestBKpivot2U(t *testing.T) {
	Ldata := [][]float64{
		[]float64{7.0, 6.0, 5.0, 9.0, 3.0, 2.0, 1.0},
		[]float64{1.0, 6.0, 9.0, 4.0, 3.0, 2.0, 1.0},
		[]float64{1.0, 2.0, 5.0, 4.0, 3.0, 5.0, 1.0},
		[]float64{1.0, 2.0, 3.0, 4.0, 3.0, 2.0, 1.0},
		[]float64{1.0, 5.0, 3.0, 4.0, 3.0, 2.0, 1.0},
		[]float64{1.0, 2.0, 3.0, 4.0, 9.0, 2.0, 1.0},
		[]float64{1.0, 2.0, 3.0, 9.0, 5.0, 6.0, 1.0}}

	Bdata := [][]float64{
		[]float64{10.0, 20.0},
		[]float64{10.0, 20.0},
		[]float64{10.0, 20.0},
		[]float64{10.0, 20.0},
		[]float64{10.0, 20.0},
		[]float64{10.0, 20.0},
		[]float64{10.0, 20.0}}

	A := matrix.FloatMatrixFromTable(Ldata, matrix.RowOrder)
	X := matrix.FloatMatrixFromTable(Bdata, matrix.RowOrder)
	N := A.Rows()
	B := matrix.FloatZeros(N, 2)
	MultSym(B, A, X, 1.0, 0.0, UPPER|LEFT)
	//t.Logf("initial B:\n%v\n", B)

	nb := 0
	W := matrix.FloatWithValue(A.Rows(), 5, 0.0)
	ipiv := make([]int, N, N)
	L, _ := DecomposeBK(A.Copy(), W, ipiv, UPPER, nb)
	t.Logf("ipiv: %v\n", ipiv)
	t.Logf("L:\n%v\n", L)

	ipiv0 := make([]int, N, N)
	nb = 4
	L0, _ := DecomposeBK(A.Copy(), W, ipiv0, UPPER, nb)
	t.Logf("ipiv: %v\n", ipiv0)
	t.Logf("L:\n%v\n", L0)
	//B0 := B.Copy()
	//SolveBK(B0, L0, ipiv0, LOWER)
	//t.Logf("B0:\n%v\n", B0)

	ipiv2 := make([]int32, N, N)
	lapack.SytrfFloat(A, ipiv2, linalg.OptUpper)
	t.Logf("ipiv2: %v\n", ipiv2)
	t.Logf("lapack A:\n%v\n", A)
	//lapack.Sytrs(A, B, ipiv2, linalg.OptUpper)
	//t.Logf("lapack B:\n%v\n", B)
	//t.Logf("B == B0: %v\n", B.AllClose(B0))
}

func _TestBKpivot2n1U(t *testing.T) {
	Ldata := [][]float64{
		[]float64{7.0, 6.0, 5.0, 9.0, 3.0, 2.0, 1.0},
		[]float64{1.0, 6.0, 9.0, 4.0, 3.0, 2.0, 1.0},
		[]float64{1.0, 2.0, 5.0, 4.0, 3.0, 5.0, 1.0},
		[]float64{1.0, 2.0, 3.0, 4.0, 3.0, 2.0, 1.0},
		[]float64{1.0, 5.0, 3.0, 4.0, 3.0, 2.0, 1.0},
		[]float64{1.0, 2.0, 3.0, 4.0, 9.0, 2.0, 1.0},
		[]float64{1.0, 2.0, 3.0, 9.0, 5.0, 6.0, 1.0}}

	Bdata := [][]float64{
		[]float64{10.0, 20.0},
		[]float64{10.0, 20.0},
		[]float64{10.0, 20.0},
		[]float64{10.0, 20.0},
		[]float64{10.0, 20.0},
		[]float64{10.0, 20.0},
		[]float64{10.0, 20.0}}

	A := matrix.FloatMatrixFromTable(Ldata, matrix.RowOrder)
	X := matrix.FloatMatrixFromTable(Bdata, matrix.RowOrder)
	N := A.Rows()
	B := matrix.FloatZeros(N, 2)
	MultSym(B, A, X, 1.0, 0.0, UPPER|LEFT)
	//t.Logf("initial B:\n%v\n", B)

	nb := 0
	W := matrix.FloatWithValue(A.Rows(), 5, 0.0)
	ipiv := make([]int, N, N)
	L, _ := DecomposeBK(A.Copy(), W, ipiv, UPPER, nb)
	t.Logf("ipiv: %v\n", ipiv)
	t.Logf("L:\n%v\n", L)

	//ipiv0 := make([]int, N, N)
	//nb = 4
	//L0, _ := DecomposeBK(A.Copy(), W, ipiv0, UPPER, nb)
	//t.Logf("ipiv: %v\n", ipiv0)
	//t.Logf("L:\n%v\n", L0)
	//B0 := B.Copy()
	//SolveBK(B0, L0, ipiv0, LOWER)
	//t.Logf("B0:\n%v\n", B0)

	ipiv2 := make([]int32, N, N)
	lapack.SytrfFloat(A, ipiv2, linalg.OptUpper)
	t.Logf("ipiv2: %v\n", ipiv2)
	t.Logf("lapack A:\n%v\n", A)
	//lapack.Sytrs(A, B, ipiv2, linalg.OptUpper)
	//t.Logf("lapack B:\n%v\n", B)
	//t.Logf("B == B0: %v\n", B.AllClose(B0))
}

func _TestBKpivot2n2U(t *testing.T) {
	Ldata := [][]float64{
		[]float64{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
		[]float64{1.0, 2.0, 2.0, 2.0, 5.0, 2.0, 2.0},
		[]float64{1.0, 2.0, 3.0, 3.0, 3.0, 3.0, 3.0},
		[]float64{1.0, 2.0, 3.0, 4.0, 4.0, 4.0, 10.0},
		[]float64{1.0, 5.0, 3.0, 4.0, 5.0, 9.0, 5.0},
		[]float64{1.0, 2.0, 3.0, 4.0, 9.0, 6.0, 6.0},
		[]float64{1.0, 2.0, 3.0, 10.0, 5.0, 6.0, 7.0}}

	Bdata := [][]float64{
		[]float64{10.0, 20.0},
		[]float64{10.0, 20.0},
		[]float64{10.0, 20.0},
		[]float64{10.0, 20.0},
		[]float64{10.0, 20.0},
		[]float64{10.0, 20.0},
		[]float64{10.0, 20.0}}

	A := matrix.FloatMatrixFromTable(Ldata, matrix.RowOrder)
	X := matrix.FloatMatrixFromTable(Bdata, matrix.RowOrder)
	N := A.Rows()
	B := matrix.FloatZeros(N, 2)
	MultSym(B, A, X, 1.0, 0.0, UPPER|LEFT)
	//t.Logf("initial B:\n%v\n", B)

	nb := 0
	W := matrix.FloatWithValue(A.Rows(), 5, 0.0)
	ipiv := make([]int, N, N)
	L, _ := DecomposeBK(A.Copy(), W, ipiv, UPPER, nb)
	t.Logf("ipiv: %v\n", ipiv)
	t.Logf("L:\n%v\n", L)

	//ipiv0 := make([]int, N, N)
	//nb = 4
	//L0, _ := DecomposeBK(A.Copy(), W, ipiv0, UPPER, nb)
	//t.Logf("ipiv: %v\n", ipiv0)
	//t.Logf("L:\n%v\n", L0)
	//B0 := B.Copy()
	//SolveBK(B0, L0, ipiv0, UPPER)
	//t.Logf("B0:\n%v\n", B0)

	ipiv2 := make([]int32, N, N)
	lapack.SytrfFloat(A, ipiv2, linalg.OptUpper)
	t.Logf("ipiv2: %v\n", ipiv2)
	t.Logf("lapack A:\n%v\n", A)
	//lapack.Sytrs(A, B, ipiv2, linalg.OptUpper)
	//t.Logf("lapack B:\n%v\n", B)
	//t.Logf("B == B0: %v\n", B.AllClose(B0))
}

func _TestBKpivot1U(t *testing.T) {
	Ldata := [][]float64{
		[]float64{7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0},
		[]float64{1.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0},
		[]float64{1.0, 2.0, 5.0, 4.0, 3.0, 5.0, 1.0},
		[]float64{1.0, 2.0, 3.0, 4.0, 3.0, 2.0, 1.0},
		[]float64{1.0, 5.0, 3.0, 4.0, 3.0, 2.0, 1.0},
		[]float64{1.0, 2.0, 3.0, 4.0, 5.0, 2.0, 1.0},
		[]float64{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 1.0}}

	Bdata := [][]float64{
		[]float64{10.0, 20.0},
		[]float64{10.0, 20.0},
		[]float64{10.0, 20.0},
		[]float64{10.0, 20.0},
		[]float64{10.0, 20.0},
		[]float64{10.0, 20.0},
		[]float64{10.0, 20.0}}

	A := matrix.FloatMatrixFromTable(Ldata, matrix.RowOrder)
	X := matrix.FloatMatrixFromTable(Bdata, matrix.RowOrder)
	N := A.Rows()
	B := matrix.FloatZeros(N, 2)
	MultSym(B, A, X, 1.0, 0.0, UPPER|LEFT)
	//t.Logf("initial B:\n%v\n", B)
	//N := 8
	//A := matrix.FloatUniformSymmetric(N)

	W := matrix.FloatWithValue(A.Rows(), 5, 0.0)

	nb := 0
	ipiv := make([]int, N, N)
	L, _ := DecomposeBK(A.Copy(), W, ipiv, UPPER, nb)
	t.Logf("ipiv: %v\n", ipiv)
	t.Logf("L:\n%v\n", L)

	nb = 4
	//ipiv0 := make([]int, N, N)
	//L0, _ := DecomposeBK(A.Copy(), W, ipiv0, LOWER, nb)
	//t.Logf("ipiv: %v\n", ipiv0)
	//t.Logf("L:\n%v\n", L0)
	//B0 := B.Copy()
	//SolveBK(B0, L0, ipiv0, LOWER)
	//t.Logf("B0:\n%v\n", B0)

	ipiv2 := make([]int32, N, N)
	lapack.SytrfFloat(A, ipiv2, linalg.OptUpper)
	t.Logf("ipiv2: %v\n", ipiv2)
	t.Logf("lapack A:\n%v\n", A)
	//lapack.Sytrs(A, B, ipiv2, linalg.OptLower)
	//t.Logf("lapack B:\n%v\n", B)
	//t.Logf("B == B0: %v\n", B.AllClose(B0))
}

func _TestBKpivot1n2U(t *testing.T) {
	Ldata := [][]float64{
		[]float64{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
		[]float64{1.0, 2.0, 2.0, 2.0, 5.0, 2.0, 2.0},
		[]float64{1.0, 2.0, 3.0, 3.0, 3.0, 3.0, 3.0},
		[]float64{1.0, 2.0, 3.0, 4.0, 4.0, 7.0, 4.0},
		[]float64{1.0, 5.0, 3.0, 4.0, 5.0, 5.0, 5.0},
		[]float64{1.0, 2.0, 3.0, 7.0, 5.0, 6.0, 6.0},
		[]float64{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0}}

	Bdata := [][]float64{
		[]float64{10.0, 20.0},
		[]float64{10.0, 20.0},
		[]float64{10.0, 20.0},
		[]float64{10.0, 20.0},
		[]float64{10.0, 20.0},
		[]float64{10.0, 20.0},
		[]float64{10.0, 20.0}}

	A := matrix.FloatMatrixFromTable(Ldata, matrix.RowOrder)
	X := matrix.FloatMatrixFromTable(Bdata, matrix.RowOrder)
	N := A.Rows()
	B := matrix.FloatZeros(N, 2)
	MultSym(B, A, X, 1.0, 0.0, UPPER|LEFT)
	t.Logf("initial B:\n%v\n", B)

	nb := 0
	W := matrix.FloatWithValue(A.Rows(), 5, 0.0)
	ipiv := make([]int, N, N)
	L, _ := DecomposeBK(A.Copy(), W, ipiv, UPPER, nb)
	t.Logf("ipiv: %v\n", ipiv)
	t.Logf("L:\n%v\n", L)

	ipiv0 := make([]int, N, N)
	nb = 4
	L0, _ := DecomposeBK(A.Copy(), W, ipiv0, UPPER, nb)
	t.Logf("ipiv0: %v\n", ipiv0)
	t.Logf("L0:\n%v\n", L0)
	//B0 := B.Copy()
	//SolveBK(B0, L0, ipiv0, LOWER)
	//t.Logf("B0:\n%v\n", B0)

	ipiv2 := make([]int32, N, N)
	lapack.Sytrf(A, ipiv2, linalg.OptUpper)
	t.Logf("ipiv2: %v\n", ipiv2)
	t.Logf("lapack A:\n%v\n", A)
	//lapack.Sytrs(A, B, ipiv2, linalg.OptLower)
	//t.Logf("lapack B:\n%v\n", B)
	//t.Logf("B == B0: %v\n", B.AllClose(B0))
}

func _TestBK2U(t *testing.T) {
	Bdata := [][]float64{
		[]float64{10.0, 20.0},
		[]float64{10.0, 20.0},
		[]float64{10.0, 20.0},
		[]float64{10.0, 20.0},
		[]float64{10.0, 20.0},
		[]float64{10.0, 20.0},
		[]float64{10.0, 20.0}}

	N := 7

	A0 := matrix.FloatNormal(N, N)
	A := matrix.FloatZeros(N, N)
	// A is symmetric, posivite definite
	Mult(A, A0, A0, 1.0, 1.0, TRANSB)

	X := matrix.FloatMatrixFromTable(Bdata, matrix.RowOrder)
	B := matrix.FloatZeros(N, 2)
	MultSym(B, A, X, 1.0, 0.0, LOWER|LEFT)
	t.Logf("initial B:\n%v\n", B)

	nb := 0
	W := matrix.FloatWithValue(A.Rows(), 5, 1.0)
	A.SetAt(4, 1, A.GetAt(4, 1)+1.0)
	A.SetAt(1, 4, A.GetAt(4, 1))

	ipiv := make([]int, N, N)
	L, _ := DecomposeBK(A.Copy(), W, ipiv, LOWER, nb)
	t.Logf("ipiv: %v\n", ipiv)
	t.Logf("L:\n%v\n", L)

	ipiv0 := make([]int, N, N)
	nb = 4
	L0, _ := DecomposeBK(A.Copy(), W, ipiv0, LOWER, nb)
	t.Logf("ipiv: %v\n", ipiv0)
	t.Logf("L:\n%v\n", L0)
	B0 := B.Copy()
	SolveBK(B0, L0, ipiv0, LOWER)
	t.Logf("B0:\n%v\n", B0)

	ipiv2 := make([]int32, N, N)
	lapack.Sytrf(A, ipiv2, linalg.OptLower)
	t.Logf("ipiv2: %v\n", ipiv2)
	t.Logf("lapack A:\n%v\n", A)
	lapack.Sytrs(A, B, ipiv2, linalg.OptLower)
	t.Logf("lapack B:\n%v\n", B)
	t.Logf("B == B0: %v\n", B.AllClose(B0))
}

func _TestBKSolveU(t *testing.T) {
	Ldata := [][]float64{
		[]float64{1.0, 2.0, 3.0, 4.0},
		[]float64{2.0, 2.0, 3.0, 4.0},
		[]float64{3.0, 3.0, 3.0, 4.0},
		[]float64{4.0, 4.0, 4.0, 4.0}}
	Xdata := [][]float64{
		[]float64{1.0, 2.0},
		[]float64{1.0, 2.0},
		[]float64{1.0, 2.0},
		[]float64{1.0, 2.0}}

	A := matrix.FloatMatrixFromTable(Ldata, matrix.RowOrder)
	X := matrix.FloatMatrixFromTable(Xdata, matrix.RowOrder)
	N := A.Rows()
	B := matrix.FloatZeros(N, 2)
	Mult(B, A, X, 1.0, 0.0, NOTRANS)
	S := matrix.FloatZeros(N, 2)
	MultSym(S, A, X, 1.0, 0.0, LOWER|UPPER)
	_, _ = S, B
	//t.Logf("B:\n%v\n", B)
	//t.Logf("S:\n%v\n", S)

	nb := 0
	W := matrix.FloatWithValue(A.Rows(), 5, 0.0)

	ipiv := make([]int, N, N)
	L, _ := DecomposeBK(A.Copy(), W, ipiv, UPPER, nb)
	t.Logf("ipiv: %v\n", ipiv)
	t.Logf("L:\n%v\n", L)
	//B0 := B.Copy()
	//SolveBK(B0, L, ipiv, LOWER)
	//t.Logf("B0:\n%v\n", B0)

	ipiv2 := make([]int32, N, N)
	lapack.Sytrf(A, ipiv2, linalg.OptUpper)
	t.Logf("ipiv2: %v\n", ipiv2)
	t.Logf("lapack A:\n%v\n", A)
	//lapack.Sytrs(A, B, ipiv2, linalg.OptLower)
	//t.Logf("lapack B:\n%v\n", B)
}

// Local Variables:
// tab-width: 4
// indent-tabs-mode: nil
// End:
