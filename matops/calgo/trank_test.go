package calgo

import (
	"github.com/henrylee2cn/algorithm/linalg"
	"github.com/henrylee2cn/algorithm/linalg/blas"
	"github.com/henrylee2cn/algorithm/matrix"
	"testing"
)

func _TestRankSmall(t *testing.T) {
	bM := 5
	bN := 5
	//bP := 5
	Adata := [][]float64{
		[]float64{1.0, 1.0, 1.0, 1.0, 1.0},
		[]float64{2.0, 2.0, 2.0, 2.0, 2.0},
		[]float64{3.0, 3.0, 3.0, 3.0, 3.0},
		[]float64{4.0, 4.0, 4.0, 4.0, 4.0},
		[]float64{5.0, 5.0, 5.0, 5.0, 5.0}}

	A := matrix.FloatMatrixFromTable(Adata, matrix.RowOrder)
	A0 := matrix.FloatMatrixFromTable(Adata, matrix.RowOrder)
	X := matrix.FloatVector([]float64{1.0, 2.0, 3.0, 4.0, 5.0})
	Y := matrix.FloatWithValue(bN, 1, 2.0)

	Ar := A.FloatArray()
	Xr := X.FloatArray()
	Yr := Y.FloatArray()

	blas.GerFloat(X, Y, A0, 1.0)

	DRankMV(Ar, Xr, Yr, 1.0, A.LeadingIndex(), 1, 1, 0, bN, 0, bM, 4, 4)
	ok := A0.AllClose(A)
	t.Logf("A0 == A1: %v\n", ok)
	if !ok {
		t.Logf("blas ger:\n%v\n", A0)
		t.Logf("A1: \n%v\n", A)
	}
}

func _TestRank(t *testing.T) {
	bM := M * 100
	bN := N * 100
	//bP := 5

	A := matrix.FloatWithValue(bM, bN, 1.0)
	A0 := matrix.FloatWithValue(bM, bN, 1.0)
	X := matrix.FloatNormal(bM, 1)
	Y := matrix.FloatNormal(bN, 1)

	Ar := A.FloatArray()
	Xr := X.FloatArray()
	Yr := Y.FloatArray()

	blas.GerFloat(X, Y, A0, 1.0)

	DRankMV(Ar, Xr, Yr, 1.0, A.LeadingIndex(), 1, 1, 0, bN, 0, bM, 4, 4)
	t.Logf("GER A0 == A1: %v\n", A0.AllClose(A))
}

func _TestMultSyrSmall(t *testing.T) {
	bN := 7
	//A := matrix.FloatNormal(bN, bN)
	//B := matrix.FloatNormal(bN, bP)
	//A := matrix.FloatWithValue(bM, bP, 1.0)
	X := matrix.FloatWithValue(bN, 1, 1.0)
	C0 := matrix.FloatZeros(bN, bN)
	C1 := matrix.FloatZeros(bN, bN)
	for i := 0; i < bN; i++ {
		X.Add(1.0+float64(i), i)
	}
	t.Logf("X=\n%v\n", X)

	Xr := X.FloatArray()
	C1r := C1.FloatArray()

	blas.SyrFloat(X, C0, 1.0, linalg.OptUpper)

	DSymmRankMV(C1r, Xr, 1.0, UPPER, C1.LeadingIndex(), 1, 0, bN, 4)
	ok := C0.AllClose(C1)
	t.Logf("UPPER C0 == C1: %v\n", ok)
	if !ok {
		t.Logf("blas: C0\n%v\n", C0)
		t.Logf("C1: C1 = A*X\n%v\n", C1)
	}

	blas.SyrFloat(X, C0, 1.0, linalg.OptLower)

	DSymmRankMV(C1r, Xr, 1.0, LOWER, C1.LeadingIndex(), 1, 0, bN, 4)
	ok = C0.AllClose(C1)
	t.Logf("LOWER C0 == C1: %v\n", ok)
	if !ok {
		t.Logf("blas: C0\n%v\n", C0)
		t.Logf("C1: C1 = A*X\n%v\n", C1)
	}
}

func _TestMultSyr2Small(t *testing.T) {
	bN := 7
	//A := matrix.FloatNormal(bN, bN)
	//B := matrix.FloatNormal(bN, bP)
	//A := matrix.FloatWithValue(bM, bP, 1.0)
	X := matrix.FloatWithValue(bN, 1, 1.0)
	Y := matrix.FloatWithValue(bN, 1, 1.0)
	C0 := matrix.FloatZeros(bN, bN)
	C1 := matrix.FloatZeros(bN, bN)
	for i := 0; i < bN; i++ {
		X.Add(1.0+float64(i), i)
		Y.Add(2.0+float64(i), i)
	}
	t.Logf("X=\n%v\nY=\n%v\n", X, Y)

	Xr := X.FloatArray()
	Yr := Y.FloatArray()
	C1r := C1.FloatArray()

	blas.Syr2Float(X, Y, C0, 1.0, linalg.OptUpper)

	DSymmRank2MV(C1r, Xr, Yr, 1.0, UPPER, C1.LeadingIndex(), 1, 1, 0, bN, 4)
	ok := C0.AllClose(C1)
	t.Logf("UPPER C0 == C1: %v\n", ok)
	if !ok {
		t.Logf("C1: C1 = A*X\n%v\n", C1)
		t.Logf("blas: C0\n%v\n", C0)
	}

	blas.Syr2Float(X, Y, C0, 1.0, linalg.OptLower)

	DSymmRank2MV(C1r, Xr, Yr, 1.0, LOWER, C1.LeadingIndex(), 1, 1, 0, bN, 4)
	ok = C0.AllClose(C1)
	t.Logf("LOWER C0 == C1: %v\n", ok)
	if !ok {
		t.Logf("blas: C0\n%v\n", C0)
		t.Logf("C1: C1 = A*X\n%v\n", C1)
	}
}

func syrkTest(t *testing.T, C, A *matrix.FloatMatrix, flags Flags, vlen, nb int) bool {
	//var B0 *matrix.FloatMatrix
	P := A.Cols()
	S := 0
	E := C.Rows()
	C0 := C.Copy()

	trans := linalg.OptNoTrans
	if flags&TRANSA != 0 {
		trans = linalg.OptTrans
		P = A.Rows()
	}
	uplo := linalg.OptUpper
	if flags&LOWER != 0 {
		uplo = linalg.OptLower
	}

	blas.SyrkFloat(A, C0, 1.0, 1.0, uplo, trans)
	if A.Rows() < 8 {
		//t.Logf("..A\n%v\n", A)
		t.Logf("  BLAS C0:\n%v\n", C0)
	}

	Ar := A.FloatArray()
	Cr := C.FloatArray()
	DSymmRankBlk(Cr, Ar, 1.0, 1.0, flags, C.LeadingIndex(), A.LeadingIndex(),
		P, S, E, vlen, nb)
	result := C0.AllClose(C)
	t.Logf("   C0 == C: %v\n", result)
	if A.Rows() < 8 {
		t.Logf("  DMRank C:\n%v\n", C)
	}
	return result
}

func _TestSyrkSmall(t *testing.T) {
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

	Adata := [][]float64{
		[]float64{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0},
		[]float64{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0}}
	A := matrix.FloatMatrixFromTable(Adata)
	t.Logf("-- SYRK UPPER --")
	syrkTest(t, U.Copy(), A, UPPER, 4, 2)
	t.Logf("-- SYRK LOWER --")
	syrkTest(t, L.Copy(), A, LOWER, 4, 2)
	t.Logf("-- SYRK UPPER, TRANSA --")
	t.Logf("A: \n%v\n", A.Transpose())
	syrkTest(t, U.Copy(), A.Transpose(), UPPER|TRANSA, 4, 2)
	t.Logf("-- SYRK LOWER --")
	syrkTest(t, L.Copy(), A.Transpose(), LOWER|TRANSA, 4, 2)
}

func syrk2Test(t *testing.T, C, A, B *matrix.FloatMatrix, flags Flags, vlen, nb int) bool {
	//var B0 *matrix.FloatMatrix
	P := A.Cols()
	S := 0
	E := C.Rows()
	C0 := C.Copy()

	trans := linalg.OptNoTrans
	if flags&TRANSA != 0 {
		trans = linalg.OptTrans
		P = A.Rows()
	}
	uplo := linalg.OptUpper
	if flags&LOWER != 0 {
		uplo = linalg.OptLower
	}

	blas.Syr2kFloat(A, B, C0, 1.0, 1.0, uplo, trans)
	if A.Rows() < 8 {
		//t.Logf("..A\n%v\n", A)
		t.Logf("  BLAS C0:\n%v\n", C0)
	}

	Ar := A.FloatArray()
	Br := B.FloatArray()
	Cr := C.FloatArray()
	DSymmRank2Blk(Cr, Ar, Br, 1.0, 1.0, flags, C.LeadingIndex(), A.LeadingIndex(),
		B.LeadingIndex(), P, S, E, vlen, nb)
	result := C0.AllClose(C)
	t.Logf("   C0 == C: %v\n", result)
	if A.Rows() < 8 {
		t.Logf("  DMRank2 C:\n%v\n", C)
	}
	return result
}

func _TestSyrk2Small(t *testing.T) {
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

	Adata := [][]float64{
		[]float64{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0},
		[]float64{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0}}
	Bdata := [][]float64{
		[]float64{7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0},
		[]float64{7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0}}
	_ = Bdata
	A := matrix.FloatMatrixFromTable(Adata)
	//B := matrix.FloatMatrixFromTable(Bdata);
	B := matrix.FloatNormal(7, 2)
	t.Logf("-- SYR2K UPPER --")
	syrk2Test(t, U.Copy(), A, B, UPPER, 4, 2)
	t.Logf("-- SYR2K LOWER --")
	syrk2Test(t, L.Copy(), A, B, LOWER, 4, 2)
	t.Logf("-- SYR2K UPPER, TRANSA --")
	//t.Logf("A: \n%v\n", A.Transpose())
	syrk2Test(t, U.Copy(), A.Transpose(), B.Transpose(), UPPER|TRANSA, 4, 2)
	t.Logf("-- SYR2K LOWER, TRANS --")
	syrk2Test(t, L.Copy(), A.Transpose(), B.Transpose(), LOWER|TRANSA, 4, 2)
}

// Local Variables:
// tab-width: 4
// indent-tabs-mode: nil
// End:
