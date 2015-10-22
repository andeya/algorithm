package calgo

import (
	"github.com/henrylee2cn/algorithm/linalg"
	"github.com/henrylee2cn/algorithm/linalg/blas"
	"github.com/henrylee2cn/algorithm/matrix"
	"testing"
)

func solveMVTest(t *testing.T, A, X0 *matrix.FloatMatrix, flags Flags, bN, bNB int) {
	X1 := X0.Copy()

	uplo := linalg.OptUpper
	diag := linalg.OptNonUnit
	if flags&LOWER != 0 {
		uplo = linalg.OptLower
	}
	if flags&UNIT != 0 {
		diag = linalg.OptUnit
	}

	blas.TrsvFloat(A, X0, uplo, diag)

	Ar := A.FloatArray()
	Xr := X1.FloatArray()
	if bN == bNB {
		DSolveUnblkMV(Xr, Ar, flags, 1, A.LeadingIndex(), bN)
	} else {
		DSolveBlkMV(Xr, Ar, flags, 1, A.LeadingIndex(), bN, bNB)
	}
	ok := X1.AllClose(X0)
	t.Logf("X1 == X0: %v\n", ok)
	if !ok && bN < 8 {
		t.Logf("A=\n%v\n", A)
		t.Logf("X0=\n%v\n", X0)
		t.Logf("blas: X0\n%v\n", X0)
		t.Logf("X1:\n%v\n", X1)
	}
}

func _TestSolveSmall(t *testing.T) {

	L := matrix.FloatMatrixFromTable(lower5, matrix.RowOrder)
	N := L.Rows()
	nb := 0

	U := L.Transpose()
	X0 := matrix.FloatWithValue(L.Rows(), 1, 1.0)
	X1 := X0.Copy()
	xsum := 0.0
	for i := 0; i < N; i++ {
		xsum += float64(i)
		X0.Add(xsum, i)
		X1.Add(xsum, -(i + 1))
	}

	t.Logf("-- SOLVE LOWER NON-UNIT ---\n")
	solveMVTest(t, L, X0.Copy(), LOWER, N, nb)
	t.Logf("-- SOLVE LOWER UNIT ---\n")
	solveMVTest(t, L, X1.Copy(), LOWER|UNIT, N, nb)

	t.Logf("-- SOLVE UPPER NON-UNIT  ---\n")
	solveMVTest(t, U, X0.Copy(), UPPER, N, nb)
	t.Logf("-- SOLVE UPPER UNIT ---\n")
	solveMVTest(t, U, X1.Copy(), UPPER|UNIT, N, nb)
}

func _TestSolveBlockedSmall(t *testing.T) {

	L := matrix.FloatMatrixFromTable(lower5, matrix.RowOrder)
	N := L.Rows()
	nb := 4

	U := L.Transpose()
	X0 := matrix.FloatWithValue(L.Rows(), 1, 1.0)
	X1 := X0.Copy()
	xsum := 0.0
	for i := 0; i < N; i++ {
		xsum += float64(i)
		X0.Add(xsum, i)
		X1.Add(xsum, -(i + 1))
	}

	t.Logf("-- SOLVE LOWER, NON-UNIT ---\n")
	solveMVTest(t, L, X0.Copy(), LOWER, N, nb)
	t.Logf("-- SOLVE UPPER, NON-UNIT ---\n")
	solveMVTest(t, U, X1.Copy(), UPPER, N, nb)

	t.Logf("-- SOLVE LOWER, UNIT ---\n")
	solveMVTest(t, L, X0.Copy(), LOWER|UNIT, N, nb)
	t.Logf("-- SOLVE UPPER, UNIT ---\n")
	solveMVTest(t, U, X1.Copy(), UPPER|UNIT, N, nb)
}

func _TestSolveRandom(t *testing.T) {
	bN := 22
	nb := 4
	L := matrix.FloatNormalSymmetric(bN, matrix.Lower)
	U := L.Transpose()
	X0 := matrix.FloatWithValue(L.Rows(), 1, 1.0)

	t.Logf("-- BLOCKED SOLVE LOWER NON-UNIT ---\n")
	solveMVTest(t, L, X0.Copy(), LOWER, bN, nb)
	t.Logf("-- BLOCKED SOLVE LOWER UNIT ---\n")
	solveMVTest(t, L, X0.Copy(), LOWER|UNIT, bN, nb)
	t.Logf("-- BLOCKED SOLVE UPPER NON-UNIT ---\n")
	solveMVTest(t, U, X0.Copy(), UPPER, bN, nb)
	t.Logf("-- BLOCKED SOLVE UPPER UNIT ---\n")
	solveMVTest(t, U, X0.Copy(), UPPER|UNIT, bN, nb)
}

func trsmSolve(t *testing.T, A *matrix.FloatMatrix, flags Flags, rand bool, nrhs, nb int) bool {
	var B0 *matrix.FloatMatrix
	side := linalg.OptLeft
	trans := linalg.OptNoTrans
	N := A.Cols()
	S := 0
	E := A.Rows()
	_ = S
	_ = E
	if flags&RIGHT != 0 {
		if rand {
			B0 = matrix.FloatNormal(nrhs, A.Rows())
		} else {
			B0 = matrix.FloatWithValue(nrhs, A.Rows(), 2.0)
		}
		side = linalg.OptRight
		E = B0.Rows()
	} else {
		if rand {
			B0 = matrix.FloatNormal(A.Rows(), nrhs)
		} else {
			B0 = matrix.FloatWithValue(A.Rows(), nrhs, 2.0)
		}
		E = B0.Cols()
	}
	B1 := B0.Copy()
	diag := linalg.OptNonUnit
	if flags&UNIT != 0 {
		diag = linalg.OptUnit
	}
	uplo := linalg.OptUpper
	if flags&LOWER != 0 {
		uplo = linalg.OptLower
	}
	if flags&TRANSA != 0 {
		trans = linalg.OptTransA
	}
	blas.TrsmFloat(A, B0, 1.0, uplo, diag, side, trans)

	Ar := A.FloatArray()
	Br := B1.FloatArray()
	if nb == 0 || nb == N {
		DSolveUnblk(Br, Ar, 1.0, flags, B1.LeadingIndex(), A.LeadingIndex(), N, S, E)
	} else {
		DSolveBlk(Br, Ar, 1.0, flags, B1.LeadingIndex(), A.LeadingIndex(), N, S, E, nb)
	}
	result := B1.AllClose(B0)
	t.Logf("B1 == B0: %v\n", result)
	if !result {
		if nrhs < 10 {
			t.Logf("blas: B0\n%v\n", B0)
			t.Logf("B1:\n%v\n", B1)
		} else {
			b0 := B0.FloatArray()
			b1 := B1.FloatArray()
			for k := 0; k < len(b0); k++ {
				if !isClose(b0[k], b1[k]) {
					t.Logf("first divergences at %d ... col %d, row %d\n", k, k/B0.Rows(), k%B0.Rows())
					break
				}
			}
		}
	}
	return result
}

func _TestTrsmSmall(t *testing.T) {
	//bN := 7
	Udata3 := [][]float64{
		[]float64{2.0, 2.0, 2.0},
		[]float64{0.0, 3.0, 3.0},
		[]float64{0.0, 0.0, 4.0}}

	Udata := [][]float64{
		[]float64{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
		[]float64{0.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0},
		[]float64{0.0, 0.0, 3.0, 3.0, 3.0, 3.0, 3.0},
		[]float64{0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 4.0},
		[]float64{0.0, 0.0, 0.0, 0.0, 5.0, 5.0, 5.0},
		[]float64{0.0, 0.0, 0.0, 0.0, 0.0, 6.0, 6.0},
		[]float64{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.0}}
	U := matrix.FloatMatrixFromTable(Udata, matrix.RowOrder)
	U3 := matrix.FloatMatrixFromTable(Udata3, matrix.RowOrder)
	_ = U
	_ = U3

	Ldata3 := [][]float64{
		[]float64{1.0, 0.0, 0.0},
		[]float64{1.0, 2.0, 0.0},
		[]float64{1.0, 2.0, 3.0}}

	Ldata := [][]float64{
		[]float64{1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		[]float64{1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		[]float64{1.0, 2.0, 3.0, 0.0, 0.0, 0.0, 0.0},
		[]float64{1.0, 2.0, 3.0, 4.0, 0.0, 0.0, 0.0},
		[]float64{1.0, 2.0, 3.0, 4.0, 5.0, 0.0, 0.0},
		[]float64{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 0.0},
		[]float64{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0}}
	L := matrix.FloatMatrixFromTable(Ldata, matrix.RowOrder)
	L3 := matrix.FloatMatrixFromTable(Ldata3, matrix.RowOrder)
	_ = L
	_ = L3

	nP := 2
	nb := 0
	t.Logf("-- TRSM-LOWER, NON-UNIT --")
	trsmSolve(t, L3, LOWER|LEFT, false, nP, nb)
	//trsmSolve(t, L, LOWER, false, nP, nb)
	//trsmSolve(t, L, LOWER, true, nP, nb)

	t.Logf("- TRSM-LOWER, UNIT --")
	trsmSolve(t, L3, LOWER|LEFT|UNIT, false, nP, nb)

	t.Logf("-- TRSM-UPPER, NON-UNIT --")
	trsmSolve(t, U3, UPPER|LEFT, false, nP, nb)
	//trsmSolve(t, U, UPPER, false, nP, nb)
	//trsmSolve(t, U, UPPER, true, nP, nb)

	t.Logf("-- TRSM-UPPER, UNIT --")
	trsmSolve(t, U3, UPPER|LEFT|UNIT, false, nP, nb)

	t.Logf("-- TRSM-UPPER, TRANS, NON_UNIT --")
	trsmSolve(t, U3, UPPER|LEFT|TRANSA, false, nP, nb)

	t.Logf("-- TRSM-UPPER, TRANS, UNIT --")
	trsmSolve(t, U3, UPPER|LEFT|TRANSA|UNIT, false, nP, nb)

	t.Logf("-- TRSM-LOWER, TRANS, NON-UNIT --")
	trsmSolve(t, L3, LOWER|LEFT|TRANSA, false, nP, nb)

	t.Logf("-- TRSM-LOWER, TRANS, UNIT --")
	trsmSolve(t, L3, LOWER|LEFT|TRANSA|UNIT, false, nP, nb)

	t.Logf("-- TRSM-UPPER, NON-UNIT, RIGHT ---")
	trsmSolve(t, U3, UPPER|RIGHT, false, nP, nb)

	t.Logf("-- TRSM-UPPER, UNIT, RIGHT ---")
	trsmSolve(t, U3, UPPER|UNIT|RIGHT, false, nP, nb)

	t.Logf("-- TRSM-LOWER, NON-UNIT, RIGHT ---")
	trsmSolve(t, L3, LOWER|RIGHT, false, nP, nb)

	t.Logf("-- TRSM-LOWER, UNIT, RIGHT ---")
	trsmSolve(t, L3, LOWER|UNIT|RIGHT, false, nP, nb)

	t.Logf("-- TRSM-UPPER, NON-UNIT, RIGHT, TRANS ---")
	trsmSolve(t, U3, UPPER|RIGHT|TRANSA, false, nP, nb)

	t.Logf("-- TRSM-UPPER, UNIT, RIGHT, TRANS ---")
	trsmSolve(t, U3, UPPER|UNIT|RIGHT|TRANSA, false, nP, nb)

	t.Logf("-- TRSM-LOWER, NON-UNIT, RIGHT, TRANS ---")
	trsmSolve(t, L3, LOWER|RIGHT|TRANSA, false, nP, nb)

	t.Logf("-- TRSM-LOWER, UNIT, RIGHT, TRANS ---")
	trsmSolve(t, L3, LOWER|UNIT|RIGHT|TRANSA, false, nP, nb)

	nP = 4
	nb = 2
	t.Logf("-- BLK TRSM-UPPER, NON-UNIT ---")
	//trsmSolve(t, U, UPPER, false, 2)
	trsmSolve(t, U, UPPER, true, nP, nb)

	t.Logf("-- BLK TRSM-UPPER, UNIT ---")
	//trsmSolve(t, U, UPPER, false, nP, nb)
	trsmSolve(t, U, UPPER|UNIT, true, nP, nb)

	t.Logf("-- BLK TRSM-LOWER, NON-UNIT ---")
	//trsmSolve(t, L, LOWER, false, nP, nb)
	trsmSolve(t, L, LOWER, true, nP, nb)

	t.Logf("-- BLK TRSM-LOWER, UNIT ---")
	//trsmSolve(t, L, LOWER, false, nP, nb)
	trsmSolve(t, L, LOWER|UNIT, true, nP, nb)

	t.Logf("-- BLK TRSM-UPPER, NON-UNIT, TRANS ---")
	//trsmSolve(t, U, UPPER|TRANSA, false, nP, nb)
	trsmSolve(t, U, UPPER|TRANSA, true, nP, nb)

	t.Logf("-- BLK TRSM-LOWER, NON-UNIT, TRANS ---")
	//trsmSolve(t, L, LOWER|TRANSA, false, nP, nb)
	trsmSolve(t, L, LOWER|TRANSA, true, nP, nb)

	t.Logf("-- BLK TRSM-UPPER, NON-UNIT, RIGHT ---")
	trsmSolve(t, U, UPPER|RIGHT, true, nP, nb)

	t.Logf("-- BLK TRSM-UPPER, UNIT, RIGHT ---")
	trsmSolve(t, U, UPPER|UNIT|RIGHT, true, nP, nb)

	t.Logf("-- BLK TRSM-UPPER, NON-UNIT, RIGHT, TRANSA ---")
	trsmSolve(t, U, UPPER|RIGHT|TRANSA, true, nP, nb)

	t.Logf("-- BLK TRSM-LOWER, NON-UNIT, RIGHT ---")
	trsmSolve(t, L, LOWER|RIGHT, true, nP, nb)

	t.Logf("-- BLK TRSM-LOWER, UNIT, RIGHT ---")
	trsmSolve(t, L, LOWER|UNIT|RIGHT, true, nP, nb)

	t.Logf("-- BLK TRSM-LOWER, NON-UNIT, RIGHT, TRANSA ---")
	trsmSolve(t, L, LOWER|RIGHT|TRANSA, true, nP, nb)
}

func TestTrsmUnblk(t *testing.T) {
	//bN := 7
	Udata3 := [][]float64{
		[]float64{2.0, 2.0, 2.0},
		[]float64{0.0, 3.0, 3.0},
		[]float64{0.0, 0.0, 4.0}}
	U3 := matrix.FloatMatrixFromTable(Udata3, matrix.RowOrder)
	_ = U3

	Ldata3 := [][]float64{
		[]float64{1.0, 0.0, 0.0},
		[]float64{1.0, 2.0, 0.0},
		[]float64{1.0, 2.0, 3.0}}
	L3 := matrix.FloatMatrixFromTable(Ldata3, matrix.RowOrder)
	_ = L3

	bN := 10
	nP := 7
	nb := 4
	L := matrix.FloatNormalSymmetric(bN, matrix.Lower)

	//t.Logf("-- TRSM-UPPER, TRANS, RIGHT, NON_UNIT --")
	//trsmSolve(t, U3, UPPER|TRANSA|RIGHT, false, 2, 0)
	//t.Logf("-- TRSM-UPPER, TRANS, RIGHT, UNIT --")
	//trsmSolve(t, U3, UPPER|TRANSA|UNIT|RIGHT, false, 2, 0)

	t.Logf("-- UNBLK TRSM-LOWER, TRANS, RIGHT, NON-UNIT --")
	trsmSolve(t, L, LOWER|TRANSA|RIGHT, false, nP, 0)
	t.Logf("-- BLK   TRSM-LOWER, TRANS, RIGHT, NON-UNIT --")
	trsmSolve(t, L, LOWER|TRANSA|RIGHT, false, nP, nb)

}

func TestTrsmBig(t *testing.T) {
	bN := 900
	nP := 40
	nb := 16
	L := matrix.FloatNormalSymmetric(bN, matrix.Lower)
	U := matrix.FloatNormalSymmetric(bN, matrix.Upper)
	_, _ = L, U

	t.Logf("-- BLK TRSM-UPPER, LEFT, NON_UNIT --")
	trsmSolve(t, U, UPPER|LEFT, false, nP, nb)
	t.Logf("-- BLK TRSM-UPPER, LEFT, UNIT --")
	trsmSolve(t, U, UPPER|LEFT|UNIT, false, nP, nb)

	t.Logf("-- BLK TRSM-LOWER, LEFT, NON_UNIT --")
	trsmSolve(t, L, LOWER|LEFT, false, nP, nb)
	t.Logf("-- BLK TRSM-LOWER, LEFT, UNIT --")
	trsmSolve(t, L, LOWER|LEFT|UNIT, false, nP, nb)

	t.Logf("-- BLK TRSM-UPPER, LEFT, TRANSA, NON_UNIT --")
	trsmSolve(t, U, UPPER|LEFT|TRANSA, false, nP, nb)
	t.Logf("-- BLK TRSM-UPPER, LEFT, TRANSA, UNIT --")
	trsmSolve(t, U, UPPER|LEFT|TRANSA|UNIT, false, nP, nb)

	t.Logf("-- BLK TRSM-LOWER, LEFT, TRANSA, NON_UNIT --")
	trsmSolve(t, L, LOWER|LEFT|TRANSA, false, nP, nb)
	t.Logf("-- BLK TRSM-LOWER, LEFT, TRANSA, UNIT --")
	trsmSolve(t, L, LOWER|LEFT|TRANSA|UNIT, false, nP, nb)

	t.Logf("-- BLK TRSM-UPPER, RIGHT, NON_UNIT --")
	trsmSolve(t, U, UPPER|RIGHT, false, nP, nb)
	t.Logf("-- BLK TRSM-UPPER, RIGHT, UNIT --")
	trsmSolve(t, U, UPPER|UNIT|RIGHT, false, nP, nb)

	t.Logf("-- BLK TRSM-LOWER, RIGHT, NON-UNIT --")
	trsmSolve(t, L, LOWER|RIGHT, false, nP, nb)
	t.Logf("-- BLK TRSM-LOWER, RIGHT, UNIT --")
	trsmSolve(t, L, LOWER|UNIT|RIGHT, false, nP, nb)

	t.Logf("-- BLK TRSM-UPPER, RIGHT, TRANSA, NON_UNIT --")
	trsmSolve(t, U, UPPER|RIGHT|TRANSA, false, nP, nb)
	t.Logf("-- BLK TRSM-UPPER, RIGHT, TRANSA, UNIT --")
	trsmSolve(t, U, UPPER|RIGHT|TRANSA|UNIT, false, nP, nb)

	nb = 0
	t.Logf("-- BLK TRSM-LOWER, RIGHT, TRANSA, NON_UNIT --")
	trsmSolve(t, L, LOWER|RIGHT|TRANSA, false, nP, nb)
	t.Logf("-- BLK TRSM-LOWER, RIGHT, TRANSA, UNIT --")
	trsmSolve(t, L, LOWER|RIGHT|TRANSA|UNIT, false, nP, nb)

}

// Local Variables:
// tab-width: 4
// indent-tabs-mode: nil
// End:
