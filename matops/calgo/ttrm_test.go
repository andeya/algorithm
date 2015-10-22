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

func trmvTest(t *testing.T, A *matrix.FloatMatrix, flags Flags, nb int) bool {
	N := A.Cols()
	//S := 0
	//E := A.Cols()
	X0 := matrix.FloatWithValue(A.Rows(), 1, 2.0)
	X1 := X0.Copy()

	trans := linalg.OptNoTrans
	if flags&TRANS != 0 {
		trans = linalg.OptTrans
	}
	diag := linalg.OptNonUnit
	if flags&UNIT != 0 {
		diag = linalg.OptUnit
	}
	uplo := linalg.OptUpper
	if flags&LOWER != 0 {
		uplo = linalg.OptLower
	}

	blas.TrmvFloat(A, X0, uplo, diag, trans)

	Ar := A.FloatArray()
	Xr := X1.FloatArray()
	if nb == 0 {
		DTrimvUnblkMV(Xr, Ar, flags, 1, A.LeadingIndex(), N)
	}
	result := X0.AllClose(X1)
	t.Logf("   X0 == X1: %v\n", result)
	if !result && A.Rows() < 8 {
		t.Logf("  BLAS TRMV X0:\n%v\n", X0)
		t.Logf("  DTrmv X1:\n%v\n", X1)
	}
	return result
}

func TestTrmvSmall(t *testing.T) {

	L := matrix.FloatMatrixFromTable(lower5, matrix.RowOrder)
	U := L.Transpose()
	nb := 0

	t.Logf("-- TRMV, LOWER --")
	trmvTest(t, L, LOWER, nb)
	t.Logf("-- TRMV, UPPER --")
	trmvTest(t, U, UPPER, nb)

	t.Logf("-- TRMV, LOWER, UNIT --")
	trmvTest(t, L, LOWER|UNIT, nb)
	t.Logf("-- TRMV, UPPER, UNIT --")
	trmvTest(t, U, UPPER|UNIT, nb)

	t.Logf("-- TRMV, LOWER, TRANS --")
	trmvTest(t, L, LOWER|TRANS, nb)
	t.Logf("-- TRMV, UPPER, TRANS --")
	trmvTest(t, U, UPPER|TRANS, nb)

	t.Logf("-- TRMV, LOWER, TRANS, UNIT --")
	trmvTest(t, L, LOWER|UNIT|TRANS, nb)
	t.Logf("-- TRMV, UPPER, TRANS, UNIT --")
	trmvTest(t, U, UPPER|UNIT|TRANS, nb)
}

func trmmTest(t *testing.T, A *matrix.FloatMatrix, flags Flags, nb int) bool {
	var B0 *matrix.FloatMatrix
	N := A.Cols()
	S := 0
	E := A.Cols()
	side := linalg.OptLeft
	if flags&RIGHT != 0 {
		B0 = matrix.FloatWithValue(2, A.Rows(), 2.0)
		side = linalg.OptRight
		E = B0.Rows()
	} else {
		B0 = matrix.FloatWithValue(A.Rows(), 2, 2.0)
		E = B0.Cols()
	}
	B1 := B0.Copy()

	trans := linalg.OptNoTrans
	if flags&TRANSA != 0 {
		trans = linalg.OptTransA
	}
	diag := linalg.OptNonUnit
	if flags&UNIT != 0 {
		diag = linalg.OptUnit
	}
	uplo := linalg.OptUpper
	if flags&LOWER != 0 {
		uplo = linalg.OptLower
	}

	blas.TrmmFloat(A, B0, 1.0, uplo, diag, trans, side)
	if A.Rows() < 8 {
		//t.Logf("..A\n%v\n", A)
		t.Logf("  BLAS B0:\n%v\n", B0)
	}

	Ar := A.FloatArray()
	Br := B1.FloatArray()
	if nb != 0 {
		DTrmmBlk(Br, Ar, 1.0, flags, B1.LeadingIndex(), A.LeadingIndex(),
			N, S, E, nb)
	} else {
		DTrmmUnblk(Br, Ar, 1.0, flags, B1.LeadingIndex(), A.LeadingIndex(),
			N, S, E, 0)
	}
	result := B0.AllClose(B1)
	t.Logf("   B0 == B1: %v\n", result)
	if A.Rows() < 8 {
		t.Logf("  DTrmm B1:\n%v\n", B1)
	}
	return result
}

func _TestTrmmUnblkSmall(t *testing.T) {

	U := matrix.FloatMatrixFromTable(upper7, matrix.RowOrder)
	U3 := matrix.FloatMatrixFromTable(upper3, matrix.RowOrder)
	_ = U
	_ = U3

	L := matrix.FloatMatrixFromTable(lower7, matrix.RowOrder)
	L3 := matrix.FloatMatrixFromTable(lower3, matrix.RowOrder)
	_ = L

	t.Logf("-- TRMM-UPPER, NON-UNIT ---")
	fail(t, trmmTest(t, U3, UPPER, 0))
	t.Logf("-- TRMM-UPPER, UNIT ---")
	fail(t, trmmTest(t, U3, UPPER|UNIT, 0))
	t.Logf("-- TRMM-UPPER, NON-UNIT, TRANSA ---")
	fail(t, trmmTest(t, U3, UPPER|TRANSA, 0))
	t.Logf("-- TRMM-UPPER, UNIT, TRANSA ---")
	fail(t, trmmTest(t, U3, UPPER|TRANSA|UNIT, 0))
	t.Logf("-- TRMM-LOWER, NON-UNIT ---")
	fail(t, trmmTest(t, L3, LOWER, 0))
	t.Logf("-- TRMM-LOWER, UNIT ---")
	fail(t, trmmTest(t, L3, LOWER|UNIT, 0))
	t.Logf("-- TRMM-LOWER, NON-UNIT, TRANSA ---")
	fail(t, trmmTest(t, L3, LOWER|TRANSA, 0))
	t.Logf("-- TRMM-LOWER, UNIT, TRANSA ---")
	fail(t, trmmTest(t, L3, LOWER|TRANSA|UNIT, 0))

	t.Logf("-- TRMM-UPPER, NON-UNIT, RIGHT ---")
	fail(t, trmmTest(t, U3, UPPER|RIGHT, 0))
	t.Logf("-- TRMM-UPPER, UNIT, RIGHT ---")
	fail(t, trmmTest(t, U3, UPPER|UNIT|RIGHT, 0))

	t.Logf("-- TRMM-LOWER, NON-UNIT, RIGHT ---")
	fail(t, trmmTest(t, L3, LOWER|RIGHT, 0))
	t.Logf("-- TRMM-LOWER, UNIT, RIGHT ---")
	fail(t, trmmTest(t, L3, LOWER|UNIT|RIGHT, 0))

	t.Logf("-- TRMM-UPPER, NON-UNIT, RIGHT, TRANSA ---")
	fail(t, trmmTest(t, U3, UPPER|RIGHT|TRANSA, 0))
	t.Logf("-- TRMM-UPPER, UNIT, RIGHT, TRANSA ---")
	fail(t, trmmTest(t, U3, UPPER|UNIT|RIGHT|TRANSA, 0))

	t.Logf("-- TRMM-LOWER, NON-UNIT, RIGHT, TRANSA ---")
	fail(t, trmmTest(t, L3, LOWER|RIGHT|TRANSA, 0))
	t.Logf("-- TRMM-LOWER, UNIT, RIGHT, TRANSA ---")
	fail(t, trmmTest(t, L3, LOWER|UNIT|RIGHT|TRANSA, 0))

}

func _TestTrmmBlkSmall(t *testing.T) {
	//bN := 7
	Udata := [][]float64{
		[]float64{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
		[]float64{0.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0},
		[]float64{0.0, 0.0, 3.0, 3.0, 3.0, 3.0, 3.0},
		[]float64{0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 4.0},
		[]float64{0.0, 0.0, 0.0, 0.0, 5.0, 5.0, 5.0},
		[]float64{0.0, 0.0, 0.0, 0.0, 0.0, 6.0, 6.0},
		[]float64{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.0}}
	U := matrix.FloatMatrixFromTable(Udata, matrix.RowOrder)
	_ = U

	Ldata := [][]float64{
		[]float64{1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		[]float64{1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		[]float64{1.0, 2.0, 3.0, 0.0, 0.0, 0.0, 0.0},
		[]float64{1.0, 2.0, 3.0, 4.0, 0.0, 0.0, 0.0},
		[]float64{1.0, 2.0, 3.0, 4.0, 5.0, 0.0, 0.0},
		[]float64{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 0.0},
		[]float64{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0}}
	L := matrix.FloatMatrixFromTable(Ldata, matrix.RowOrder)
	_ = L

	t.Logf("-- TRMM-UPPER, NON-UNIT ---")
	fail(t, trmmTest(t, U, UPPER, 2))
	t.Logf("-- TRMM-UPPER, NON-UNIT, TRANSA ---")
	fail(t, trmmTest(t, U, UPPER|TRANSA, 2))
	t.Logf("-- TRMM-LOWER, NON-UNIT ---")
	fail(t, trmmTest(t, L, LOWER, 2))
	t.Logf("-- TRMM-LOWER, NON-UNIT, TRANSA ---")
	fail(t, trmmTest(t, L, LOWER|TRANSA, 2))
	t.Logf("-- TRMM-UPPER, RIGHT, NON-UNIT ---")
	fail(t, trmmTest(t, U, UPPER|RIGHT, 2))
	t.Logf("-- TRMM-UPPER, RIGHT, NON-UNIT, TRANSA ---")
	fail(t, trmmTest(t, U, UPPER|RIGHT|TRANSA, 2))
	t.Logf("-- TRMM-LOWER, RIGHT, NON-UNIT ---")
	fail(t, trmmTest(t, U, LOWER|RIGHT, 2))
	t.Logf("-- TRMM-LOWER, RIGHT, NON-UNIT, TRANSA ---")
	fail(t, trmmTest(t, U, LOWER|RIGHT|TRANSA, 2))
}

// Local Variables:
// tab-width: 4
// indent-tabs-mode: nil
// End:
