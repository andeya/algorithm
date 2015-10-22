// Copyright (c) Harri Rautila, 2012

// This file is part of go.opt package. It is free software, distributed
// under the terms of GNU Lesser General Public License Version 3, or any later
// version. See the COPYING tile included in this archive.

package cvx

import (
	"github.com/henrylee2cn/algorithm/cvx/sets"
	"github.com/henrylee2cn/algorithm/matrix"
	"testing"
)

func TestConeQp(t *testing.T) {
	adata := [][]float64{
		[]float64{0.3, -0.4, -0.2, -0.4, 1.3},
		[]float64{0.6, 1.2, -1.7, 0.3, -0.3},
		[]float64{-0.3, 0.0, 0.6, -1.2, -2.0}}

	// reference values from cvxopt coneqp.py
	xref := []float64{0.72558318685981904, 0.61806264311119252, 0.30253527966423444}
	sref := []float64{
		0.72558318685981904, 0.61806264311119263,
		0.30253527966423449, 1.00000000000041678,
		-0.72558318686012169, -0.61806264311145032,
		-0.30253527966436067}
	zref := []float64{
		0.00000003332583626, 0.00000005116586239,
		0.00000009993673262, 0.56869648433154019,
		0.41264857754144563, 0.35149286573190930,
		0.17201618570052318}

	A := matrix.FloatMatrixFromTable(adata, matrix.ColumnOrder)
	b := matrix.FloatVector([]float64{1.5, 0.0, -1.2, -0.7, 0.0})

	_, n := A.Size()
	N := n + 1 + n

	h := matrix.FloatZeros(N, 1)
	h.SetIndex(n, 1.0)

	I0 := matrix.FloatDiagonal(n, -1.0)
	I1 := matrix.FloatIdentity(n)
	G, _ := matrix.FloatMatrixStacked(matrix.StackDown, I0, matrix.FloatZeros(1, n), I1)

	At := A.Transpose()
	P := matrix.Times(At, A)
	q := matrix.Times(At, b).Scale(-1.0)

	dims := sets.DSetNew("l", "q", "s")
	dims.Set("l", []int{n})
	dims.Set("q", []int{n + 1})

	var solopts SolverOptions
	solopts.MaxIter = 10
	solopts.ShowProgress = false
	sol, err := ConeQp(P, q, G, h, nil, nil, dims, &solopts, nil)
	if err == nil {
		fail := false
		x := sol.Result.At("x")[0]
		s := sol.Result.At("s")[0]
		z := sol.Result.At("z")[0]
		t.Logf("Optimal\n")
		t.Logf("x=\n%v\n", x.ToString("%.9f"))
		t.Logf("s=\n%v\n", s.ToString("%.9f"))
		t.Logf("z=\n%v\n", z.ToString("%.9f"))
		xe, _ := nrmError(matrix.FloatVector(xref), x)
		if xe > TOL {
			t.Logf("x differs [%.3e] from exepted too much.", xe)
			fail = true
		}
		se, _ := nrmError(matrix.FloatVector(sref), s)
		if se > TOL {
			t.Logf("s differs [%.3e] from exepted too much.", se)
			fail = true
		}
		ze, _ := nrmError(matrix.FloatVector(zref), z)
		if ze > TOL {
			t.Logf("z differs [%.3e] from exepted too much.", ze)
			fail = true
		}
		if fail {
			t.Fail()
		}
	}

}

// Local Variables:
// tab-width: 4
// End:
