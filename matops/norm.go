// Copyright (c) Harri Rautila, 2013

// This file is part of github.com/hrautila/matops package. It is free software,
// distributed under the terms of GNU Lesser General Public License Version 3, or
// any later version. See the COPYING tile included in this archive.

package matops

import (
	"github.com/henrylee2cn/algorithm/matrix"
)

type Norms int

const (
	NORM_ONE = Norms(1)
	NORM_TWO = Norms(2)
	NORM_INF = Norms(-1)
)

func mNorm1(A *matrix.FloatMatrix) float64 {
	var amax float64 = 0.0
	var col matrix.FloatMatrix
	for k := 0; k < A.Cols(); k++ {
		col.SubMatrixOf(A, 0, k, A.Rows(), 1)
		cmax := ASum(&col)
		if cmax > amax {
			amax = cmax
		}
	}
	return amax
}

func mNormInf(A *matrix.FloatMatrix) float64 {
	var amax float64 = 0.0
	var row matrix.FloatMatrix
	for k := 0; k < A.Rows(); k++ {
		row.SubMatrixOf(A, k, 0, A.Cols(), 1)
		rmax := ASum(&row)
		if rmax > amax {
			amax = rmax
		}
	}
	return amax
}

/*
 * Compute matrix and vector norms.
 *
 * Arguments
 *  X    A real valued matrix or vector
 *
 *  norm Norm to compute
 *         NORM_ONE, NORM_TWO, NORM_INF
 *
 * Note: matrix NORM_TWO not yet implemented.
 */
func NormP(X *matrix.FloatMatrix, norm Norms) float64 {
	if isVector(X) {
		switch norm {
		case NORM_ONE:
			return ASum(X)
		case NORM_TWO:
			return Norm2(X)
		case NORM_INF:
			return AMax(X)
		}
		return 0.0
	}
	switch norm {
	case NORM_ONE:
		return mNorm1(X)
	case NORM_TWO:
		return 0.0
	case NORM_INF:
		return mNormInf(X)
	}
	return 0.0
}

// Local Variables:
// tab-width: 4
// indent-tabs-mode: nil
// End:
