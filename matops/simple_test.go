// Copyright (c) Harri Rautila, 2012,2013

// This file is part of github.com/hrautila/matops package. It is free software,
// distributed under the terms of GNU Lesser General Public License Version 3, or
// any later version. See the COPYING tile included in this archive.

package matops

import (
	"github.com/henrylee2cn/algorithm/matrix"
	"testing"
)

func TestCompile(t *testing.T) {
	t.Logf("Compile OK\n")
}

func newPivots(sz int) *pPivots {
	p := new(pPivots)
	p.pivots = make([]int, sz, sz)
	return p
}

func (p *pPivots) Size() int {
	return len(p.pivots)
}

func _TestSwap(t *testing.T) {
	A := matrix.FloatUniform(5, 5)
	t.Logf("A:\n%v\n", A)
	swapRows(A, 0, 4)
	t.Logf("A: 0,4 swapped\n%v\n", A)
	swapRows(A, 4, 0)
	t.Logf("A: 4,0 swapped again\n%v\n", A)
}

func _TestPivot1Dv(t *testing.T) {
	var pT, pB, p0, p1, p2 pPivots
	p := newPivots(6)
	partitionPivot2x1(&pT, &pB, p, 0, pBOTTOM)
	t.Logf("m(pT)=%d, m(pB)=%d\n", pT.Size(), pB.Size())
	k := 0
	for pT.Size() < p.Size() {
		t.Logf("m(pB)=%d; %v\n", pB.Size(), pB)
		repartPivot2x1to3x1(&pT, &p0, &p1, &p2, p, 1, pBOTTOM)
		p1.pivots[0] += k
		t.Logf("m(p0)=%d, m(p2)=%d, p1=%d\n", p0.Size(), p2.Size(), p1.pivots[0])
		contPivot3x1to2x1(&pT, &pB, &p0, &p1, p, pBOTTOM)
		k += 1
	}
	t.Logf("pivots: %v\n", p.pivots)
}

func _TestPartition1Dh(t *testing.T) {
	var AL, AR, A0, a1, A2 matrix.FloatMatrix
	A := matrix.FloatZeros(1, 6)
	partition1x2(&AL, &AR, A, 0, pRIGHT)
	t.Logf("m(AL)=%d, m(AR)=%d\n", AL.Cols(), AR.Cols())
	for AL.Cols() < A.Cols() {
		AR.Add(1.0)
		t.Logf("m(AR)=%d; %v\n", AR.Cols(), AR)
		repartition1x2to1x3(&AL, &A0, &a1, &A2, A, 1, pRIGHT)
		t.Logf("m(A0)=%d, m(A2)=%d, a1=%.1f\n", A0.Cols(), A2.Cols(), a1.Float())
		continue1x3to1x2(&AL, &AR, &A0, &a1, A, pRIGHT)
	}
}

func _TestPartition2D(t *testing.T) {
	var ATL, ATR, ABL, ABR, As matrix.FloatMatrix
	var A00, a01, A02, a10, a11, a12, A20, a21, A22 matrix.FloatMatrix

	A := matrix.FloatZeros(6, 6)
	As.SubMatrixOf(A, 1, 1, 4, 4)
	As.SetIndexes(1.0)
	partition2x2(&ATL, &ATR, &ABL, &ABR, &As, 0)
	t.Logf("ATL:\n%v\n", &ATL)

	for ATL.Rows() < As.Rows() {
		repartition2x2to3x3(&ATL,
			&A00, &a01, &A02,
			&a10, &a11, &a12,
			&A20, &a21, &A22, &As, 1)
		t.Logf("m(a12)=%d [%d], m(a11)=%d\n", a12.Cols(), a12.NumElements(), a11.NumElements())
		a11.Add(1.0)
		a21.Add(-2.0)

		continue3x3to2x2(&ATL, &ATR, &ABL, &ABR, &A00, &a11, &A22, &As)
	}
	t.Logf("A:\n%v\n", A)
}

func _TestViewUpdate(t *testing.T) {
	Adata2 := [][]float64{
		[]float64{4.0, 2.0, 2.0},
		[]float64{6.0, 4.0, 2.0},
		[]float64{4.0, 6.0, 1.0},
	}

	A := matrix.FloatMatrixFromTable(Adata2, matrix.RowOrder)
	N := A.Rows()

	// simple LU decomposition without pivoting
	var A11, a10, a01, a00 matrix.FloatMatrix
	for k := 1; k < N; k++ {
		a00.SubMatrixOf(A, k-1, k-1, 1, 1)
		a01.SubMatrixOf(A, k-1, k, 1, A.Cols()-k)
		a10.SubMatrixOf(A, k, k-1, A.Rows()-k, 1)
		A11.SubMatrixOf(A, k, k)
		//t.Logf("A11: %v  a01: %v\n", A11, a01)
		a10.Scale(1.0 / a00.Float())
		MVRankUpdate(&A11, &a10, &a01, -1.0)
	}

	Ld := TriLU(A.Copy())
	Ud := TriU(A)
	t.Logf("Ld:\n%v\nUd:\n%v\n", Ld, Ud)
	An := matrix.FloatZeros(N, N)
	Mult(An, Ld, Ud, 1.0, 1.0, NOTRANS)
	t.Logf("A == Ld*Ud: %v\n", An.AllClose(An))
}

// Local Variables:
// tab-width: 4
// indent-tabs-mode: nil
// End:
