// Copyright (c) Harri Rautila, 2012,2013

// This file is part of github.com/hrautila/matops package. It is free software,
// distributed under the terms of GNU Lesser General Public License Version 3, or
// any later version. See the COPYING tile included in this archive.

package matops

import (
	"github.com/henrylee2cn/algorithm/matrix"
	"testing"
	//"math"
)

func TestUpperCHOL(t *testing.T) {
	N := 60
	K := 30
	nb := 16
	Z := matrix.FloatUniform(N, N)
	A := matrix.Times(Z, Z.Transpose())
	B := matrix.FloatUniform(N, K)
	X := B.Copy()

	// R = chol(A) = U.T*U
	R, _ := DecomposeCHOL(TriU(A.Copy()), UPPER, nb)

	// X = A.-1*B = U.-1*(U.-T*B)
	SolveCHOL(X, R, UPPER)

	// B = B - A*X
	Mult(B, A, X, -1.0, 1.0, NONE)

	// ||B - A*X||_1
	nrm := NormP(B, NORM_ONE)
	t.Logf("||B - A*X||_1: %e\n", nrm)
}

func TestLowerCHOL(t *testing.T) {
	N := 60
	K := 30
	nb := 16
	Z := matrix.FloatUniform(N, N)
	A := matrix.Times(Z, Z.Transpose())
	B := matrix.FloatUniform(N, K)
	X := B.Copy()

	// R = chol(A) = L*L.T
	R, _ := DecomposeCHOL(A.Copy(), LOWER, nb)

	// X = A.-1*B = L.-T*(L.-1*B)
	SolveCHOL(X, R, LOWER)

	// B = B - A*X
	Mult(B, A, X, -1.0, 1.0, NONE)

	nrm := NormP(B, NORM_ONE)
	t.Logf("||B - A*X||_1: %e\n", nrm)
}

// Local Variables:
// tab-width: 4
// indent-tabs-mode: nil
// End:
