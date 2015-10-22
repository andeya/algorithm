// Copyright (c) Harri Rautila, 2012,2013

// This file is part of github.com/hrautila/matops package. It is free software,
// distributed under the terms of GNU Lesser General Public License Version 3, or
// any later version. See the COPYING tile included in this archive.

package matops

import (
	"github.com/henrylee2cn/algorithm/matrix"
	"testing"
)

func TestLU(t *testing.T) {
	N := 60
	K := 30
	nb := 12
	A := matrix.FloatUniform(N, N)
	B := matrix.FloatUniform(N, K)
	X := B.Copy()
	piv := make([]int, N, N)

	// R = lu(A) = P*L*U
	R, _ := DecomposeLU(A.Copy(), piv, nb)

	// X = A.-1*B = U.-1*(L.-1*B)
	SolveLU(X, R, piv, NONE)

	// B = B - A*X
	Mult(B, A, X, -1.0, 1.0, NONE)

	nrm := NormP(B, NORM_ONE)
	t.Logf("||B - A*X||_1: %e\n", nrm)
}

// Local Variables:
// tab-width: 4
// indent-tabs-mode: nil
// End:
