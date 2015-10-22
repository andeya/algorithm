// Copyright (c) Harri Rautila, 2013

// This file is part of github.com/hrautila/matops package. It is free software,
// distributed under the terms of GNU Lesser General Public License Version 3, or
// any later version. See the COPYING tile included in this archive.

package matops

import (
	"errors"
	"github.com/henrylee2cn/algorithm/matrix"
	//"math"
	//"fmt"
)

func imin(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func m(A *matrix.FloatMatrix) int {
	return A.Rows()
}

func n(A *matrix.FloatMatrix) int {
	return A.Cols()
}

var decompNB int = 0

// Set global decomposition block size for blocked versions.
func DecomposeBlockSize(nb int) {
	decompNB = nb
}

// unblocked LU decomposition w/o pivots, FLAME LU nopivots variant 5
func unblockedLUnoPiv(A *matrix.FloatMatrix) (err error) {
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

		// a21 = a21/a11
		//a21.Scale(1.0/a11.Float())
		InvScale(&a21, a11.Float())
		// A22 = A22 - a21*a12
		err = MVRankUpdate(&A22, &a21, &a12, -1.0)

		continue3x3to2x2(
			&ATL, &ATR,
			&ABL, &ABR, &A00, &a11, &A22, A, pBOTTOMRIGHT)
	}
	return
}

// blocked LU decomposition w/o pivots, FLAME LU nopivots variant 5
func blockedLUnoPiv(A *matrix.FloatMatrix, nb int) (err error) {
	var ATL, ATR, ABL, ABR matrix.FloatMatrix
	var A00, A01, A02, A10, A11, A12, A20, A21, A22 matrix.FloatMatrix

	err = nil
	partition2x2(
		&ATL, &ATR,
		&ABL, &ABR, A, 0, 0, pTOPLEFT)

	for ATL.Rows() < A.Rows() {
		repartition2x2to3x3(&ATL,
			&A00, &A01, &A02,
			&A10, &A11, &A12,
			&A20, &A21, &A22, A, nb, pBOTTOMRIGHT)

		// A00 = LU(A00)
		unblockedLUnoPiv(&A11)
		// A12 = trilu(A00)*A12.-1  (TRSM)
		SolveTrm(&A12, &A11, 1.0, LEFT|LOWER|UNIT)
		// A21 = A21.-1*triu(A00) (TRSM)
		SolveTrm(&A21, &A11, 1.0, RIGHT|UPPER)
		// A22 = A22 - A21*A12
		Mult(&A22, &A21, &A12, -1.0, 1.0, NOTRANS)

		continue3x3to2x2(
			&ATL, &ATR,
			&ABL, &ABR, &A00, &A11, &A22, A, pBOTTOMRIGHT)
	}
	return
}

// unblocked LU decomposition with pivots: FLAME LU variant 3
func unblockedLUpiv(A *matrix.FloatMatrix, p *pPivots) error {
	var err error
	var ATL, ATR, ABL, ABR matrix.FloatMatrix
	var A00, a01, A02, a10, a11, a12, A20, a21, A22 matrix.FloatMatrix
	var AL, AR, A0, a1, A2, aB1, AB0 matrix.FloatMatrix
	var pT, pB, p0, p1, p2 pPivots

	err = nil
	partition2x2(
		&ATL, &ATR,
		&ABL, &ABR, A, 0, 0, pTOPLEFT)
	partition1x2(
		&AL, &AR, A, 0, pLEFT)
	partitionPivot2x1(
		&pT,
		&pB, p, 0, pTOP)

	for ATL.Rows() < A.Rows() && ATL.Cols() < A.Cols() {
		repartition2x2to3x3(&ATL,
			&A00, &a01, &A02,
			&a10, &a11, &a12,
			&A20, &a21, &A22 /**/, A, 1, pBOTTOMRIGHT)
		repartition1x2to1x3(&AL,
			&A0, &a1, &A2 /**/, A, 1, pRIGHT)
		repartPivot2x1to3x1(&pT,
			&p0, &p1, &p2 /**/, p, 1, pBOTTOM)

		// apply previously computed pivots
		applyPivots(&a1, &p0)

		// a01 = trilu(A00) \ a01 (TRSV)
		MVSolveTrm(&a01, &A00, 1.0, LOWER|UNIT)
		// a11 = a11 - a10 *a01
		a11.Add(Dot(&a10, &a01, -1.0))
		// a21 = a21 -A20*a01
		MVMult(&a21, &A20, &a01, -1.0, 1.0, NOTRANS)

		// pivot index on current column [a11, a21].T
		aB1.SubMatrixOf(&ABR, 0, 0, ABR.Rows(), 1)
		pivotIndex(&aB1, &p1)

		// pivots to current column
		applyPivots(&aB1, &p1)

		// a21 = a21 / a11
		InvScale(&a21, a11.Float())

		// apply pivots to previous columns
		AB0.SubMatrixOf(&ABL, 0, 0)
		applyPivots(&AB0, &p1)
		// scale last pivots to origin matrix row numbers
		p1.pivots[0] += ATL.Rows()

		continue3x3to2x2(
			&ATL, &ATR,
			&ABL, &ABR, &A00, &a11, &A22, A, pBOTTOMRIGHT)
		continue1x3to1x2(
			&AL, &AR, &A0, &a1, A, pRIGHT)
		contPivot3x1to2x1(
			&pT,
			&pB, &p0, &p1, p, pBOTTOM)
	}
	if ATL.Cols() < A.Cols() {
		applyPivots(&ATR, p)
		SolveTrm(&ATR, &ATL, 1.0, LEFT|UNIT|LOWER)
	}
	return err
}

// blocked LU decomposition with pivots: FLAME LU variant 3
func blockedLUpiv(A *matrix.FloatMatrix, p *pPivots, nb int) error {
	var err error
	var ATL, ATR, ABL, ABR matrix.FloatMatrix
	var A00, A01, A02, A10, A11, A12, A20, A21, A22 matrix.FloatMatrix
	var AL, AR, A0, A1, A2, AB1, AB0 matrix.FloatMatrix
	var pT, pB, p0, p1, p2 pPivots

	err = nil
	partition2x2(
		&ATL, &ATR,
		&ABL, &ABR, A, 0, 0, pTOPLEFT)
	partition1x2(
		&AL, &AR, A, 0, pLEFT)
	partitionPivot2x1(
		&pT,
		&pB, p, 0, pTOP)

	for ATL.Rows() < A.Rows() && ATL.Cols() < A.Cols() {
		repartition2x2to3x3(&ATL,
			&A00, &A01, &A02,
			&A10, &A11, &A12,
			&A20, &A21, &A22, A, nb, pBOTTOMRIGHT)
		repartition1x2to1x3(&AL,
			&A0, &A1, &A2 /**/, A, nb, pRIGHT)
		repartPivot2x1to3x1(&pT,
			&p0, &p1, &p2 /**/, p, nb, pBOTTOM)

		// apply previously computed pivots
		applyPivots(&A1, &p0)

		// a01 = trilu(A00) \ a01 (TRSV)
		SolveTrm(&A01, &A00, 1.0, LOWER|UNIT)
		// A11 = A11 - A10*A01
		Mult(&A11, &A10, &A01, -1.0, 1.0, NOTRANS)
		// A21 = A21 - A20*A01
		Mult(&A21, &A20, &A01, -1.0, 1.0, NOTRANS)

		// LU_piv(AB1, p1)
		AB1.SubMatrixOf(&ABR, 0, 0, ABR.Rows(), A11.Cols())
		unblockedLUpiv(&AB1, &p1)

		// apply pivots to previous columns
		AB0.SubMatrixOf(&ABL, 0, 0)
		applyPivots(&AB0, &p1)
		// scale last pivots to origin matrix row numbers
		for k, _ := range p1.pivots {
			p1.pivots[k] += ATL.Rows()
		}

		continue3x3to2x2(
			&ATL, &ATR,
			&ABL, &ABR /**/, &A00, &A11, &A22, A, pBOTTOMRIGHT)
		continue1x3to1x2(
			&AL, &AR /**/, &A0, &A1, A, pRIGHT)
		contPivot3x1to2x1(
			&pT,
			&pB /**/, &p0, &p1, p, pBOTTOM)
	}
	if ATL.Cols() < A.Cols() {
		applyPivots(&ATR, p)
		SolveTrm(&ATR, &ATL, 1.0, LEFT|UNIT|LOWER)
	}
	return err
}

/*
 * Compute an LU factorization of a general M-by-N matrix using
 * partial pivoting with row interchanges.
 *
 * Arguments:
 *   A      On entry, the M-by-N matrix to be factored. On exit the factors
 *          L and U from factorization A = P*L*U, the unit diagonal elements
 *          of L are not stored.
 *
 *   pivots On exit the pivot indices.
 *
 *   nb     Blocking factor for blocked invocations. If bn == 0 or
 *          min(M,N) < nb unblocked algorithm is used.
 *
 * Returns:
 *  LU factorization and error indicator.
 *
 * Compatible with lapack.DGETRF
 */
func DecomposeLU(A *matrix.FloatMatrix, pivots []int, nb int) (*matrix.FloatMatrix, error) {
	var err error
	mlen := imin(A.Rows(), A.Cols())
	if len(pivots) < mlen {
		return A, errors.New("pivot array < min(A.Rows(),A.Cols())")
	}
	// clear pivot array
	for k, _ := range pivots {
		pivots[k] = 0
	}
	if mlen <= nb || nb == 0 {
		err = unblockedLUpiv(A, &pPivots{pivots})
	} else {
		err = blockedLUpiv(A, &pPivots{pivots}, nb)
	}
	return A, err
}

/*
 * Compute an LU factorization of a general M-by-N matrix without pivoting.
 *
 * Arguments:
 *   A   On entry, the M-by-N matrix to be factored. On exit the factors
 *       L and U from factorization A = P*L*U, the unit diagonal elements
 *       of L are not stored.
 *
 *   nb  Blocking factor for blocked invocations. If bn == 0 or
 *       min(M,N) < nb unblocked algorithm is used.
 *
 * Returns:
 *  LU factorization and error indicator.
 *
 * Compatible with lapack.DGETRF
 */
func DecomposeLUnoPiv(A *matrix.FloatMatrix, nb int) (*matrix.FloatMatrix, error) {
	var err error
	mlen := imin(A.Rows(), A.Cols())
	if mlen <= nb || nb == 0 {
		err = unblockedLUnoPiv(A)
	} else {
		err = blockedLUnoPiv(A, nb)
	}
	return A, err
}

/*
 * Solve a system of linear equations A*X = B or A.T*X = B with general M-by-N
 * matrix A using the LU factorizatoin computed by DecomposeLU().
 *
 * Arguments:
 *  B      On entry, the right hand side matrix B. On exit, the solution matrix X.
 *
 *  A      The factor L and U from the factorization A = P*L*U as computed by
 *         DecomposeLU()
 *
 *  pivots The pivot indices from DecomposeLU().
 *
 *  flags  The indicator of the form of the system of equations.
 *         If flags&TRANSA then system is transposed. All other values
 *         indicate non transposed system.
 *
 * Compatible with lapack.DGETRS.
 */
func SolveLU(B, A *matrix.FloatMatrix, pivots []int, flags Flags) error {
	var err error = nil
	applyPivots(B, &pPivots{pivots})
	if flags&TRANSA != 0 {
		// transposed X = A.-1*B == (L.T*U.T).-1*B == U.-T*(L.-T*B)
		SolveTrm(B, A, 1.0, LOWER|UNIT|TRANSA)
		SolveTrm(B, A, 1.0, UPPER|TRANSA)
	} else {
		// non-transposed X = A.-1*B == (L*U).-1*B == U.-1*(L.-1*B)
		SolveTrm(B, A, 1.0, LOWER|UNIT)
		SolveTrm(B, A, 1.0, UPPER)
	}

	return err
}

// Local Variables:
// tab-width: 4
// indent-tabs-mode: nil
// End:
