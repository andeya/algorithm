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
 *  ( a11  a12 )   ( 1   0   )( d1  0   )( l  l21.t )
 *  ( a21  A22 )   ( l21 L22 )(  0  A22 )( 0  L22.t )
 *
 *   a11  =   d1
 *   a21  =   l21*d1                       => l21 = a21/d1
 *   A22  =   l21*d1*l21.t + L22*D2*L22.t  => L22 = A22 - l21*d1*l21t
 */
func unblkLowerLDLnoPiv(A *matrix.FloatMatrix) (err error) {
	var ATL, ATR, ABL, ABR matrix.FloatMatrix
	var A00, a10, a11, A20, a21, A22 matrix.FloatMatrix

	err = nil
	partition2x2(
		&ATL, &ATR,
		&ABL, &ABR, A, 0, 0, pTOPLEFT)

	for ATL.Rows() < A.Rows() {
		repartition2x2to3x3(&ATL,
			&A00, nil, nil,
			&a10, &a11, nil,
			&A20, &a21, &A22, A, 1, pBOTTOMRIGHT)

		// --------------------------------------------------------

		// d11 = a11; no-op

		// A22 = A22 - l21*d11*l21.T = A22 - a21*a21.T/a11; triangular update
		err = MVUpdateTrm(&A22, &a21, &a21, -1.0/a11.Float(), LOWER)

		// l21 = a21/a11
		InvScale(&a21, a11.Float())
		// ---------------------------------------------------------

		continue3x3to2x2(
			&ATL, &ATR,
			&ABL, &ABR, &A00, &a11, &A22, A, pBOTTOMRIGHT)
	}
	return
}

/*
 *  ( A11  a12 )   ( U11 u12 )( D1  0  )( U11.t 0 )
 *  ( a21  a22 )   (  0   1  )(  0  d2 )( u12.t 1 )
 *
 *   a22  =   d2
 *   a01  =   u12*d2                       => u12 = a12/d2
 *   A11  =   u12*d2*u12.t + U11*D1*U11.t  => U11 = A11 - u12*d2*u12.t
 */
func unblkUpperLDLnoPiv(A *matrix.FloatMatrix) (err error) {
	var ATL, ATR, ABL, ABR matrix.FloatMatrix
	var A00, a01, A02, a11, a12, A22 matrix.FloatMatrix

	err = nil
	partition2x2(
		&ATL, &ATR,
		&ABL, &ABR, A, 0, 0, pBOTTOMRIGHT)

	for ATL.Rows() > 0 {
		repartition2x2to3x3(&ATL,
			&A00, &a01, &A02,
			nil, &a11, &a12,
			nil, nil, &A22, A, 1, pTOPLEFT)

		// --------------------------------------------------------

		// A00 = A00 - u01*d11*u01.T = A00 - a01*a01.T/a11; triangular update
		err = MVUpdateTrm(&A00, &a01, &a01, -1.0/a11.Float(), UPPER)

		// u01 = a01/a11
		InvScale(&a01, a11.Float())
		// ---------------------------------------------------------

		continue3x3to2x2(
			&ATL, &ATR,
			&ABL, &ABR, &A00, &a11, &A22, A, pTOPLEFT)
	}
	return
}

/*
 *  ( A11  A12 )   ( L11   0  )( D1  0  )( L11.t  L21.t )
 *  ( A21  A22 )   ( L21  L22 )(  0  D2 )(   0    L22.t )
 *
 *   A11  =   L11*D1*L11.t                 -> L11\D1 = LDL(A11)
 *   A12  =   L11*D1*L21.t
 *   A21  =   L21*D1*L11.t                 => L21 = A21*(D1*L11.t).-1 = A21*L11.-T*D1.-1
 *   A22  =   L21*D1*L21.t + L22*D2*L22.t  => L22 = A22 - L21*D1*L21.t
 */
func blkLowerLDLnoPiv(A, W *matrix.FloatMatrix, nb int) (err error) {
	var ATL, ATR, ABL, ABR matrix.FloatMatrix
	var A00, A10, A11, A20, A21, A22 matrix.FloatMatrix
	var D1, wrk matrix.FloatMatrix

	err = nil
	partition2x2(
		&ATL, &ATR,
		&ABL, &ABR, A, 0, 0, pTOPLEFT)

	for ATL.Rows() < A.Rows() {
		repartition2x2to3x3(&ATL,
			&A00, nil, nil,
			&A10, &A11, nil,
			&A20, &A21, &A22, A, nb, pBOTTOMRIGHT)

		// --------------------------------------------------------

		// A11 = LDL(A11)
		unblkLowerLDLnoPiv(&A11)
		A11.Diag(&D1)

		// A21 = A21*A11.-T
		SolveTrm(&A21, &A11, 1.0, LOWER|UNIT|RIGHT|TRANSA)
		// A21 = A21*D1.-1
		SolveDiag(&A21, &D1, RIGHT)

		// W = D1*L21.T = L21*D1
		W.SubMatrix(&wrk, 0, 0, A21.Rows(), nb)
		A21.CopyTo(&wrk)
		MultDiag(&wrk, &D1, RIGHT)

		// A22 = A22 - L21*D1*L21.T = A22 - L21*W
		UpdateTrm(&A22, &A21, &wrk, -1.0, 1.0, LOWER|TRANSB)

		// ---------------------------------------------------------

		continue3x3to2x2(
			&ATL, &ATR,
			&ABL, &ABR, &A00, &A11, &A22, A, pBOTTOMRIGHT)
	}
	return
}

func blkUpperLDLnoPiv(A, W *matrix.FloatMatrix, nb int) (err error) {
	var ATL, ATR, ABL, ABR matrix.FloatMatrix
	var A00, A01, A02, A11, A12, A22 matrix.FloatMatrix
	var D1, wrk matrix.FloatMatrix

	err = nil
	partition2x2(
		&ATL, &ATR,
		&ABL, &ABR, A, 0, 0, pBOTTOMRIGHT)

	for ATL.Rows() > 0 {
		repartition2x2to3x3(&ATL,
			&A00, &A01, &A02,
			nil, &A11, &A12,
			nil, nil, &A22, A, nb, pTOPLEFT)

		// --------------------------------------------------------

		// A11 = LDL(A11)
		unblkUpperLDLnoPiv(&A11)
		A11.Diag(&D1)

		// A01 = A01*A11.-T
		SolveTrm(&A01, &A11, 1.0, UPPER|UNIT|RIGHT|TRANSA)
		// A01 = A01*D1.-1
		SolveDiag(&A01, &D1, RIGHT)

		// W = D1*U01.T = U01*D1
		W.SubMatrix(&wrk, 0, 0, A01.Rows(), nb)
		A01.CopyTo(&wrk)
		MultDiag(&wrk, &D1, RIGHT)

		// A00 = A00 - U01*D1*U01.T = A22 - U01*W.T
		UpdateTrm(&A00, &A01, &wrk, -1.0, 1.0, UPPER|TRANSB)

		// ---------------------------------------------------------

		continue3x3to2x2(
			&ATL, &ATR,
			&ABL, &ABR, &A00, &A11, &A22, A, pTOPLEFT)
	}
	return
}

/*
 * Compute an LDLT factorization of a symmetric N-by-N matrix without pivoting.
 *
 * Arguments:
 *   A      On entry, the N-by-N matrix to be factored. On exit the factor
 *          L and 1-by-1 diagonal D from factorization A = L*D*L.T, the unit diagonal
 *          of L are not stored.
 *
 *   W      Work space for blocking invocations, matrix of size N-by-nb.
 *
 *   flags  Indicator bits.
 *
 *   nb     Blocking factor for blocked invocations. If bn == 0 or
 *          N < nb unblocked algorithm is used.
 *
 * Returns:
 *  LDL factorization and error indicator.
 *
 */
func DecomposeLDLnoPiv(A, W *matrix.FloatMatrix, flags Flags, nb int) (*matrix.FloatMatrix, error) {
	var err error
	if A.Cols() != A.Rows() {
		return nil, errors.New("A not a square matrix")
	}
	if A.Cols() < nb || nb == 0 {
		if flags&LOWER != 0 {
			err = unblkLowerLDLnoPiv(A)
		} else {
			err = unblkUpperLDLnoPiv(A)
		}
	} else {
		if W == nil {
			return nil, errors.New("No workspace for blocking invocation")
		}
		if flags&LOWER != 0 {
			err = blkLowerLDLnoPiv(A, W, nb)
		} else {
			err = blkUpperLDLnoPiv(A, W, nb)
		}
	}
	return A, err
}

/*
 * Solves a system system of linear equations A*X = B with symmetric positive
 * definite matrix A using the LDL factorization A = U.T*D*U or A = L*D*L.T
 * computed by DecomposeLDLnoPiv().
 *
 * Arguments:
 *  B   On entry, the right hand side matrix B. On exit, the solution
 *      matrix X.
 *
 *  A   The triangular factor U or L from LDL factorization as computed by
 *      DecomposeLDLnoPiv().
 *
 *  flags Indicator of which factor is stored in A. If flags&UPPER then upper
 *        triangle of A is stored. If flags&LOWER then lower triangle of A is
 *        stored.
 */
func SolveLDLnoPiv(B, A *matrix.FloatMatrix, flags Flags) {
	if flags&UPPER != 0 {
		// X = (U*D*U.T).-1*B => U.-T*(D.-1*(U.-1*B))
		SolveTrm(B, A, 1.0, UPPER|UNIT)
		SolveDiag(B, A, LEFT)
		SolveTrm(B, A, 1.0, UPPER|UNIT|TRANSA)
	} else if flags&LOWER != 0 {
		// X = (L*D*L.T).-1*B = L.-T*(D*-1(L.-1*B))
		SolveTrm(B, A, 1.0, LOWER|UNIT)
		SolveDiag(B, A, LEFT)
		SolveTrm(B, A, 1.0, LOWER|UNIT|TRANSA)
	}
}

// Local Variables:
// tab-width: 4
// indent-tabs-mode: nil
// End:
