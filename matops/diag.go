// Copyright (c) Harri Rautila, 2012,2013

// This file is part of github.com/hrautila/matops package. It is free software,
// distributed under the terms of GNU Lesser General Public License Version 3, or
// any later version. See the COPYING tile included in this archive.

package matops

import (
	"errors"
	"github.com/henrylee2cn/algorithm/matrix"
)

/*
 * Compute
 *   C = C*diag(D)      flags & RIGHT == true
 *   C = diag(D)*C      flags & LEFT  == true
 *
 * Arguments
 *   C     M-by-N matrix if flags&RIGHT == true or N-by-M matrix if flags&LEFT == true
 *
 *   D     N element column or row vector or N-by-N matrix
 *
 *   flags Indicator bits, LEFT or RIGHT
 */
func MultDiag(C, D *matrix.FloatMatrix, flags Flags) {
	var c, d0 matrix.FloatMatrix
	if D.Cols() == 1 {
		// diagonal is column vector
		switch flags & (LEFT | RIGHT) {
		case LEFT:
			// scale rows; for each column element-wise multiply with D-vector
			for k := 0; k < C.Cols(); k++ {
				C.SubMatrix(&c, 0, k, C.Rows(), 1)
				c.Mul(D)
			}
		case RIGHT:
			// scale columns
			for k := 0; k < C.Cols(); k++ {
				C.SubMatrix(&c, 0, k, C.Rows(), 1)
				// scale the column
				c.Scale(D.GetAt(k, 0))
			}
		}
	} else {
		// diagonal is row vector
		var d *matrix.FloatMatrix
		if D.Rows() == 1 {
			d = D
		} else {
			D.SubMatrix(&d0, 0, 0, 1, D.Cols(), D.LeadingIndex()+1)
			d = &d0
		}
		switch flags & (LEFT | RIGHT) {
		case LEFT:
			for k := 0; k < C.Rows(); k++ {
				C.SubMatrix(&c, k, 0, 1, C.Cols())
				// scale the row
				c.Scale(d.GetAt(0, k))
			}
		case RIGHT:
			// scale columns
			for k := 0; k < C.Cols(); k++ {
				C.SubMatrix(&c, 0, k, C.Rows(), 1)
				// scale the column
				c.Scale(d.GetAt(0, k))
			}
		}
	}
}

/*
 * Compute
 *   X = B*diag(D).-1      flags & RIGHT == true
 *   X = diag(D).-1*C      flags & LEFT  == true
 *
 * Arguments:
 *   B     M-by-N matrix if flags&RIGHT == true or N-by-M matrix if flags&LEFT == true
 *
 *   D     N element column or row vector or N-by-N matrix
 *
 *   flags Indicator bits, LEFT or RIGHT
 */
func SolveDiag(B, D *matrix.FloatMatrix, flags Flags) {
	var c, d0 matrix.FloatMatrix
	if D.Cols() == 1 {
		// diagonal is column vector
		switch flags & (LEFT | RIGHT) {
		case LEFT:
			// scale rows; for each column element-wise multiply with D-vector
			for k := 0; k < B.Cols(); k++ {
				B.SubMatrix(&c, 0, k, B.Rows(), 1)
				c.Div(D)
			}
		case RIGHT:
			// scale columns
			for k := 0; k < B.Cols(); k++ {
				B.SubMatrix(&c, 0, k, B.Rows(), 1)
				// scale the column
				c.Scale(1.0 / D.GetAt(k, 0))
			}
		}
	} else {
		var d *matrix.FloatMatrix
		if D.Rows() == 1 {
			d = D
		} else {
			D.SubMatrix(&d0, 0, 0, 1, D.Cols(), D.LeadingIndex()+1)
			d = &d0
		}
		switch flags & (LEFT | RIGHT) {
		case LEFT:
			for k := 0; k < B.Rows(); k++ {
				B.SubMatrix(&c, k, 0, 1, B.Cols())
				// scale the row
				c.Scale(1.0 / d.GetAt(0, k))
			}
		case RIGHT:
			// scale columns
			for k := 0; k < B.Cols(); k++ {
				B.SubMatrix(&c, 0, k, B.Rows(), 1)
				// scale the column
				c.Scale(1.0 / d.GetAt(0, k))
			}
		}
	}
}

/*
 * Generic rank update of diagonal matrix.
 *   diag(D) = diag(D) + alpha * x * y.T
 *
 * Arguments:
 *   D     N element column or row vector or N-by-N matrix
 *
 *   x, y  N element vectors
 *
 *   alpha scalar
 */
func MVUpdateDiag(D, x, y *matrix.FloatMatrix, alpha float64) error {
	var d *matrix.FloatMatrix
	var dvec matrix.FloatMatrix

	if !isVector(x) || !isVector(y) {
		return errors.New("x, y not vectors")
	}
	if D.Rows() > 0 && D.Cols() == D.Rows() {
		D.Diag(&dvec)
		d = &dvec
	} else if isVector(D) {
		d = D
	} else {
		return errors.New("D not a diagonal")
	}

	N := d.NumElements()
	for k := 0; k < N; k++ {
		val := d.GetIndex(k)
		val += x.GetIndex(k) * y.GetIndex(k) * alpha
		d.SetIndex(k, val)
	}
	return nil
}

// Local Variables:
// tab-width: 4
// indent-tabs-mode: nil
// End:
