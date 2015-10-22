// Copyright (c) Harri Rautila, 2013

// This file is part of github.com/hrautila/matops package. It is free software,
// distributed under the terms of GNU Lesser General Public License Version 3, or
// any later version. See the COPYING tile included in this archive.

package matops

import (
	"errors"
	"fmt"
	"github.com/henrylee2cn/algorithm/matrix"
	"math"
)

func unblockedCHOL(A *matrix.FloatMatrix, flags Flags, nr int) (err error) {
	var ATL, ATR, ABL, ABR matrix.FloatMatrix
	var A00, a01, A02, a10, a11, a12, A20, a21, A22 matrix.FloatMatrix

	err = nil
	partition2x2(
		&ATL, &ATR,
		&ABL, &ABR, A, 0, 0, pTOPLEFT)

	for ATL.Rows() < A.Rows() {
		repartition2x2to3x3(&ATL,
			&A00, &a01, &A02,
			&a10, &a11, &a12,
			&A20, &a21, &A22, A, 1, pBOTTOMRIGHT)

		// a11 = sqrt(a11)
		aval := math.Sqrt(a11.Float())
		if math.IsNaN(aval) {
			panic(fmt.Sprintf("illegal value at %d: %e", nr+ATL.Rows(), a11.Float()))
		}
		a11.SetAt(0, 0, aval)

		if flags&LOWER != 0 {
			// a21 = a21/a11
			InvScale(&a21, a11.Float())
			// A22 = A22 - a21*a21' (SYR)
			err = MVRankUpdateSym(&A22, &a21, -1.0, flags)
		} else {
			// a21 = a12/a11
			InvScale(&a12, a11.Float())
			// A22 = A22 - a12'*a12 (SYR)
			err = MVRankUpdateSym(&A22, &a12, -1.0, flags)
		}

		continue3x3to2x2(
			&ATL, &ATR,
			&ABL, &ABR, &A00, &a11, &A22, A, pBOTTOMRIGHT)
	}
	return
}

func blockedCHOL(A *matrix.FloatMatrix, flags Flags, nb int) error {
	var err error
	var ATL, ATR, ABL, ABR matrix.FloatMatrix
	var A00, A01, A02, A10, A11, A12, A20, A21, A22 matrix.FloatMatrix

	err = nil
	partition2x2(
		&ATL, &ATR,
		&ABL, &ABR, A, 0, 0, pTOPLEFT)

	for ATL.Rows() < A.Rows() && ATL.Cols() < A.Cols() {
		repartition2x2to3x3(&ATL,
			&A00, &A01, &A02,
			&A10, &A11, &A12,
			&A20, &A21, &A22, A, nb, pBOTTOMRIGHT)

		// A11 = chol(A11)
		err = unblockedCHOL(&A11, flags, ATL.Rows())

		if flags&LOWER != 0 {
			// A21 = A21 * tril(A11).-1
			SolveTrm(&A21, &A11, 1.0, RIGHT|LOWER|TRANSA)
			// A22 = A22 - A21*A21.T
			RankUpdateSym(&A22, &A21, -1.0, 1.0, LOWER)
		} else {
			// A12 = triu(A11).-1 * A12
			SolveTrm(&A12, &A11, 1.0, UPPER|TRANSA)
			// A22 = A22 - A12.T*A12
			RankUpdateSym(&A22, &A12, -1.0, 1.0, UPPER|TRANSA)
		}

		continue3x3to2x2(
			&ATL, &ATR,
			&ABL, &ABR, &A00, &A11, &A22, A, pBOTTOMRIGHT)
	}
	return err
}

/*
 * Compute the Cholesky factorization of a symmetric positive definite
 * N-by-N matrix A.
 *
 * Arguments:
 *  A     On entry, the symmetric matrix A. If flags&UPPER the upper triangular part
 *        of A contains the upper triangular part of the matrix A, and strictly
 *        lower part A is not referenced. If flags&LOWER the lower triangular part
 *        of a contains the lower triangular part of the matrix A. Likewise, the
 *        strictly upper part of A is not referenced. On exit, factor U or L from the
 *        Cholesky factorization A = U.T*U or A = L*L.T
 *
 *  flags The matrix structure indicator, UPPER for upper tridiagonal and LOWER for
 *        lower tridiagonal matrix.
 *
 *  nb    The blocking factor for blocked invocations. If nb == 0 or N < nb unblocked
 *        algorithm is used.
 *
 * Compatible with lapack.DPOTRF
 */
func DecomposeCHOL(A *matrix.FloatMatrix, flags Flags, nb int) (*matrix.FloatMatrix, error) {
	var err error
	if A.Cols() != A.Rows() {
		return A, errors.New("A not a square matrix")
	}
	if A.Cols() < nb || nb == 0 {
		err = unblockedCHOL(A, flags, 0)
	} else {
		err = blockedCHOL(A, flags, nb)
	}
	return A, err
}

/*
 * Solves a system system of linear equations A*X = B with symmetric positive
 * definite matrix A using the Cholesky factorization A = U.T*U or A = L*L.T
 * computed by DecomposeCHOL().
 *
 * Arguments:
 *  B   On entry, the right hand side matrix B. On exit, the solution
 *      matrix X.
 *
 *  A   The triangular factor U or L from Cholesky factorization as computed by
 *      DecomposeCHOL().
 *
 *  flags Indicator of which factor is stored in A. If flags&UPPER then upper
 *        triangle of A is stored. If flags&LOWER then lower triangle of A is
 *        stored.
 *
 * Compatible with lapack.DPOTRS.
 */
func SolveCHOL(B, A *matrix.FloatMatrix, flags Flags) {
	// A*X = B; X = A.-1*B == (LU).-1*B == U.-1*L.-1*B == U.-1*(L.-1*B)
	if flags&UPPER != 0 {
		// X = (U.T*U).-1*B => U.-1*(U.-T*B)
		SolveTrm(B, A, 1.0, UPPER|TRANSA)
		SolveTrm(B, A, 1.0, UPPER)
	} else if flags&LOWER != 0 {
		// X = (L*L.T).-1*B = L.-T*(L.1*B)
		SolveTrm(B, A, 1.0, LOWER)
		SolveTrm(B, A, 1.0, LOWER|TRANSA)
	}
}

// Local Variables:
// tab-width: 4
// indent-tabs-mode: nil
// End:
