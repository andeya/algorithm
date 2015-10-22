// Copyright (c) Harri Rautila, 2012,2013

// This file is part of github.com/hrautila/matops package. It is free software,
// distributed under the terms of GNU Lesser General Public License Version 3, or
// any later version. See the COPYING tile included in this archive.

package matops

import (
	"github.com/henrylee2cn/algorithm/matrix"
	"testing"
)

func TestInvUpper(t *testing.T) {
	N := 36
	nb := 12
	A := matrix.FloatUniformSymmetric(N, matrix.Upper)
	I := matrix.FloatDiagonal(N, 1.0)
	I0 := matrix.FloatZeros(N, N)

	A0 := A.Copy()
	// A0 = A.-1
	InverseTrm(A0, UPPER, 0)
	Mult(I0, A, A0, 1.0, 0.0, NOTRANS)
	I0.Minus(I)
	nrm := NormP(I0, NORM_ONE)
	t.Logf("unblk: ||I - A*A.-1||_1 : %e\n", nrm)

	A0 = A.Copy()
	// A0 = A.-1
	InverseTrm(A0, UPPER, nb)
	Mult(I0, A, A0, 1.0, 0.0, NOTRANS)

	I0.Minus(I)
	nrm = NormP(I0, NORM_ONE)
	t.Logf("blk:   ||I - A*A.-1||_1 : %e\n", nrm)
}

func TestInvLower(t *testing.T) {
	N := 36
	nb := 12
	A := matrix.FloatUniformSymmetric(N, matrix.Lower)
	I := matrix.FloatDiagonal(N, 1.0)
	I0 := matrix.FloatZeros(N, N)

	A0 := A.Copy()
	// A0 = A.-1
	InverseTrm(A0, LOWER, 0)
	Mult(I0, A, A0, 1.0, 0.0, NOTRANS)
	I0.Minus(I)
	nrm := NormP(I0, NORM_ONE)
	t.Logf("unblk: ||I - A*A.-1||_1 : %e\n", nrm)

	A0 = A.Copy()
	// A0 = A.-1
	InverseTrm(A0, LOWER, nb)
	Mult(I0, A, A0, 1.0, 0.0, NOTRANS)
	I0.Minus(I)
	nrm = NormP(I0, NORM_ONE)
	t.Logf("blk:   ||I - A*A.-1||_1 : %e\n", nrm)
}

// Local Variables:
// tab-width: 4
// indent-tabs-mode: nil
// End:
