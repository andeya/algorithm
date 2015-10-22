// Copyright (c) Harri Rautila, 2012,2013

// This file is part of github.com/hrautila/matops package. It is free software,
// distributed under the terms of GNU Lesser General Public License Version 3, or
// any later version. See the COPYING tile included in this archive.

package matops

import (
	"github.com/henrylee2cn/algorithm/matrix"
	//"fmt"
)

// Inverse NON-UNIT diagonal tridiagonal matrix
func unblockedInverseLower(A *matrix.FloatMatrix) (err error) {
	var ATL, ATR, ABL, ABR matrix.FloatMatrix
	var A00, a10t, a11, A20, a21, A22 matrix.FloatMatrix

	err = nil
	partition2x2(
		&ATL, &ATR,
		&ABL, &ABR, A, 0, 0, pTOPLEFT)

	for ATL.Rows() < A.Rows() {
		repartition2x2to3x3(&ATL,
			&A00, nil, nil,
			&a10t, &a11, nil,
			&A20, &a21, &A22, A, 1, pBOTTOMRIGHT)
		// -------------------------------------------------
		aval := a11.Float()

		// a21 = -a21/a11
		InvScale(&a21, -aval)
		// A20 = A20 + a21*a10.t
		MVRankUpdate(&A20, &a21, &a10t, 1.0)
		// a10 = a10/a11
		InvScale(&a10t, aval)
		// a11 = 1.0/a11
		a11.SetAt(0, 0, 1.0/aval)

		// -------------------------------------------------
		continue3x3to2x2(
			&ATL, &ATR,
			&ABL, &ABR, &A00, &a11, &A22, A, pBOTTOMRIGHT)
	}
	return
}

// Inverse NON-UNIT diagonal tridiagonal matrix
func unblockedInverseUpper(A *matrix.FloatMatrix) (err error) {
	var ATL, ATR, ABL, ABR matrix.FloatMatrix
	var A00, a01, A02, a11, a12t, A22 matrix.FloatMatrix

	err = nil
	partition2x2(
		&ATL, &ATR,
		&ABL, &ABR, A, 0, 0, pTOPLEFT)

	for ATL.Rows() < A.Rows() {
		repartition2x2to3x3(&ATL,
			&A00, &a01, &A02,
			nil, &a11, &a12t,
			nil, nil, &A22, A, 1, pBOTTOMRIGHT)
		// -------------------------------------------------
		aval := a11.Float()

		// a12 = -a12/a11
		InvScale(&a12t, -aval)
		// A02 = A02 + a01*a12
		MVRankUpdate(&A02, &a01, &a12t, 1.0)
		// a01 = a01/a11
		InvScale(&a01, aval)
		// a11 = 1.0/a11
		a11.SetAt(0, 0, 1.0/aval)

		// -------------------------------------------------
		continue3x3to2x2(
			&ATL, &ATR,
			&ABL, &ABR, &A00, &a11, &A22, A, pBOTTOMRIGHT)
	}
	return
}

// Inverse UNIT diagonal tridiagonal matrix
func unblockedInverseUnitLower(A *matrix.FloatMatrix) (err error) {
	var ATL, ATR, ABL, ABR matrix.FloatMatrix
	var A00, a10t, a11, A20, a21, A22 matrix.FloatMatrix

	err = nil
	partition2x2(
		&ATL, &ATR,
		&ABL, &ABR, A, 0, 0, pTOPLEFT)

	for ATL.Rows() < A.Rows() {
		repartition2x2to3x3(&ATL,
			&A00, nil, nil,
			&a10t, &a11, nil,
			&A20, &a21, &A22, A, 1, pBOTTOMRIGHT)
		// -------------------------------------------------

		// a21 = -a21
		Scale(&a21, -1.0)
		// A20 = A20 + a21*a10.t
		MVRankUpdate(&A20, &a21, &a10t, 1.0)

		// -------------------------------------------------
		continue3x3to2x2(
			&ATL, &ATR,
			&ABL, &ABR, &A00, &a11, &A22, A, pBOTTOMRIGHT)
	}
	return
}

// Inverse NON-UNIT diagonal tridiagonal matrix
func unblockedInverseUnitUpper(A *matrix.FloatMatrix) (err error) {
	var ATL, ATR, ABL, ABR matrix.FloatMatrix
	var A00, a01, A02, a11, a12t, A22 matrix.FloatMatrix

	err = nil
	partition2x2(
		&ATL, &ATR,
		&ABL, &ABR, A, 0, 0, pTOPLEFT)

	for ATL.Rows() < A.Rows() {
		repartition2x2to3x3(&ATL,
			&A00, &a01, &A02,
			nil, &a11, &a12t,
			nil, nil, &A22, A, 1, pBOTTOMRIGHT)
		// -------------------------------------------------

		// a12 = -a12/a11
		Scale(&a12t, -1.0)
		// A02 = A02 + a01*a12.t
		MVRankUpdate(&A02, &a01, &a12t, 1.0)

		// -------------------------------------------------
		continue3x3to2x2(
			&ATL, &ATR,
			&ABL, &ABR, &A00, &a11, &A22, A, pBOTTOMRIGHT)
	}
	return
}

// Inverse tridiagonal matrix
func blockedInverseLower(A *matrix.FloatMatrix, flags Flags, nb int) (err error) {
	var ATL, ATR, ABL, ABR matrix.FloatMatrix
	var A00, A10, A11, A20, A21, A22 matrix.FloatMatrix

	err = nil
	partition2x2(
		&ATL, &ATR,
		&ABL, &ABR, A, 0, 0, pTOPLEFT)

	for ATL.Rows() < A.Rows() {
		repartition2x2to3x3(&ATL,
			&A00, nil, nil,
			&A10, &A11, nil,
			&A20, &A21, &A22, A, nb, pBOTTOMRIGHT)
		// -------------------------------------------------
		// libflame, variant 3

		// A21 = -A21 * triu(A11).-1
		SolveTrm(&A21, &A11, -1.0, flags|RIGHT)
		// A210 = A20 + A21*A10
		Mult(&A20, &A21, &A10, 1.0, 1.0, NONE)
		// A10 = tri(A11).-1*A10
		SolveTrm(&A10, &A11, 1.0, flags)
		// A11 = inv(A11)
		if flags&UNIT != 0 {
			unblockedInverseUnitLower(&A11)
		} else {
			unblockedInverseLower(&A11)
		}
		// -------------------------------------------------
		continue3x3to2x2(
			&ATL, &ATR,
			&ABL, &ABR, &A00, &A11, &A22, A, pBOTTOMRIGHT)
	}
	return
}

func blockedInverseUpper(A *matrix.FloatMatrix, flags Flags, nb int) (err error) {
	var ATL, ATR, ABL, ABR matrix.FloatMatrix
	var A00, A01, A02, A11, A12, A22 matrix.FloatMatrix

	err = nil
	partition2x2(
		&ATL, &ATR,
		&ABL, &ABR, A, 0, 0, pTOPLEFT)

	for ATL.Rows() < A.Rows() {
		repartition2x2to3x3(&ATL,
			&A00, &A01, &A02,
			nil, &A11, &A12,
			nil, nil, &A22, A, nb, pBOTTOMRIGHT)
		// -------------------------------------------------
		// libflame, variant 1

		// A01 = A00*A01
		MultTrm(&A01, &A00, 1.0, flags)
		// A01 = -A01 / triu(A11)
		SolveTrm(&A01, &A11, -1.0, flags|RIGHT)
		// A11 = inv(A11)
		if flags&UNIT != 0 {
			unblockedInverseUnitUpper(&A11)
		} else {
			unblockedInverseUpper(&A11)
		}
		// -------------------------------------------------
		continue3x3to2x2(
			&ATL, &ATR,
			&ABL, &ABR, &A00, &A11, &A22, A, pBOTTOMRIGHT)
	}
	return
}

func InverseTrm(A *matrix.FloatMatrix, flags Flags, nb int) (*matrix.FloatMatrix, error) {
	var err error = nil
	if nb == 0 || A.Cols() < nb {
		if flags&UNIT != 0 {
			if flags&LOWER != 0 {
				err = unblockedInverseUnitLower(A)
			} else {
				err = unblockedInverseUnitUpper(A)
			}
		} else {
			if flags&LOWER != 0 {
				err = unblockedInverseLower(A)
			} else {
				err = unblockedInverseUpper(A)
			}
		}
	} else {
		if flags&LOWER != 0 {
			err = blockedInverseLower(A, flags, nb)
		} else {
			err = blockedInverseUpper(A, flags, nb)
		}
	}
	return A, err
}

// Local Variables:
// tab-width: 4
// indent-tabs-mode: nil
// End:
