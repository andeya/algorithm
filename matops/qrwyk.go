// Copyright (c) Harri Rautila, 2013

// This file is part of github.com/hrautila/matops package. It is free software,
// distributed under the terms of GNU Lesser General Public License Version 3, or
// any later version. See the COPYING tile included in this archive.

package matops

import (
	"errors"
	"github.com/henrylee2cn/algorithm/matrix"
	//"fmt"
)

// Build Q in place by applying elementary reflectors in reverse order to
// an implied identity matrix.  This forms Q = H(1)H(2) ... H(k)
//
// this is compatibe with lapack.DORG2R
func unblockedBuildQ(A, tau, w *matrix.FloatMatrix, kb int) error {
	var err error = nil
	var ATL, ATR, ABL, ABR matrix.FloatMatrix
	var A00, a01, A02, a10t, a11, a12t, A20, a21, A22 matrix.FloatMatrix
	var tT, tB matrix.FloatMatrix
	var t0, tau1, t2, w1 matrix.FloatMatrix
	var mb int
	var rowvec bool

	mb = A.Rows() - A.Cols()
	rowvec = tau.Rows() == 1

	partition2x2(
		&ATL, &ATR,
		&ABL, &ABR, A, mb, 0, pBOTTOMRIGHT)

	if rowvec {
		partition1x2(
			&tT, &tB, tau, 0, pRIGHT)
	} else {
		partition2x1(
			&tT,
			&tB, tau, 0, pBOTTOM)
	}

	// clearing of the columns of the right and setting ABR to unit diagonal
	// (only if not applying all reflectors, kb > 0)

	for ATL.Rows() > 0 && ATL.Cols() > 0 {
		repartition2x2to3x3(&ATL,
			&A00, &a01, &A02,
			&a10t, &a11, &a12t,
			&A20, &a21, &A22, A, 1, pTOPLEFT)
		if rowvec {
			repartition1x2to1x3(&tT,
				&t0, &tau1, &t2, tau, 1, pLEFT)
		} else {
			repartition2x1to3x1(&tT,
				&t0,
				&tau1,
				&t2, tau, 1, pTOP)
		}

		// --------------------------------------------------------

		// adjust workspace to correct size
		w.SubMatrix(&w1, 0, 0, 1, a12t.Cols())
		// apply Householder reflection from left
		applyHHTo2x1(&tau1, &a21, &a12t, &A22, &w1, LEFT)

		// apply (in-place) current elementary reflector to unit vector
		a21.Scale(-tau1.Float())
		a11.SetAt(0, 0, 1.0-tau1.Float())

		// zero the upper part
		a01.SetIndexes(0.0)

		// --------------------------------------------------------
		continue3x3to2x2(
			&ATL, &ATR,
			&ABL, &ABR, &A00, &a11, &A22, A, pTOPLEFT)
		if rowvec {
			continue1x3to1x2(
				&tT, &tB, &t0, &tau1, tau, pLEFT)
		} else {
			continue3x1to2x1(
				&tT,
				&tB, &t0, &tau1, tau, pTOP)
		}
	}
	return err
}

func blockedBuildQ(A, tau, W *matrix.FloatMatrix, nb int) error {
	var err error = nil
	var ATL, ATR, ABL, ABR, AL matrix.FloatMatrix
	var A00, A01, A02, A10, A11, A12, A20, A21, A22 matrix.FloatMatrix
	var tT, tB matrix.FloatMatrix
	var t0, tau1, t2, Tw, Wrk matrix.FloatMatrix
	var mb int

	mb = A.Rows() - A.Cols()
	Twork := matrix.FloatZeros(nb, nb)

	partition2x2(
		&ATL, &ATR,
		&ABL, &ABR, A, mb, 0, pBOTTOMRIGHT)
	partition2x1(
		&tT,
		&tB, tau, 0, pBOTTOM)

	// clearing of the columns of the right and setting ABR to unit diagonal
	// (only if not applying all reflectors, kb > 0)

	for ATL.Rows() > 0 && ATL.Cols() > 0 {
		repartition2x2to3x3(&ATL,
			&A00, &A01, &A02,
			&A10, &A11, &A12,
			&A20, &A21, &A22, A, nb, pTOPLEFT)
		repartition2x1to3x1(&tT,
			&t0,
			&tau1,
			&t2, tau, nb, pTOP)

		// --------------------------------------------------------

		// build block reflector from current block
		merge2x1(&AL, &A11, &A21)
		Twork.SubMatrix(&Tw, 0, 0, A11.Cols(), A11.Cols())
		unblkQRBlockReflector(&Tw, &AL, &tau1)

		// update with current block reflector (I - Y*T*Y.T)*Atrailing
		W.SubMatrix(&Wrk, 0, 0, A12.Cols(), A11.Cols())
		updateWithQT(&A12, &A22, &A11, &A21, &Tw, &Wrk, nb, false)

		// use unblocked version to compute current block
		W.SubMatrix(&Wrk, 0, 0, 1, A11.Cols())
		unblockedBuildQ(&AL, &tau1, &Wrk, 0)

		// zero upper part
		A01.SetIndexes(0.0)

		// --------------------------------------------------------
		continue3x3to2x2(
			&ATL, &ATR,
			&ABL, &ABR, &A00, &A11, &A22, A, pTOPLEFT)
		continue3x1to2x1(
			&tT,
			&tB, &t0, &tau1, tau, pTOP)
	}
	return err
}

func blockedBuildQT(A, T, W *matrix.FloatMatrix, nb int) error {
	var err error = nil
	var ATL, ATR, ABL, ABR, AL matrix.FloatMatrix
	var A00, A01, A11, A12, A21, A22 matrix.FloatMatrix
	var TTL, TTR, TBL, TBR matrix.FloatMatrix
	var T00, T01, T02, T11, T12, T22 matrix.FloatMatrix
	var tau1, Wrk matrix.FloatMatrix
	var mb int

	mb = A.Rows() - A.Cols()

	partition2x2(
		&ATL, &ATR,
		&ABL, &ABR, A, mb, 0, pBOTTOMRIGHT)
	partition2x2(
		&TTL, &TTR,
		&TBL, &TBR, T, 0, 0, pBOTTOMRIGHT)

	// clearing of the columns of the right and setting ABR to unit diagonal
	// (only if not applying all reflectors, kb > 0)

	for ATL.Rows() > 0 && ATL.Cols() > 0 {
		repartition2x2to3x3(&ATL,
			&A00, &A01, nil,
			nil, &A11, &A12,
			nil, &A21, &A22, A, nb, pTOPLEFT)
		repartition2x2to3x3(&TTL,
			&T00, &T01, &T02,
			nil, &T11, &T12,
			nil, nil, &T22, T, nb, pTOPLEFT)

		// --------------------------------------------------------

		// update with current block reflector (I - Y*T*Y.T)*Atrailing
		W.SubMatrix(&Wrk, 0, 0, A12.Cols(), A11.Cols())
		updateWithQT(&A12, &A22, &A11, &A21, &T11, &Wrk, nb, false)

		// use unblocked version to compute current block
		W.SubMatrix(&Wrk, 0, 0, 1, A11.Cols())
		// elementary scalar coefficients on the diagonal, column vector
		T11.Diag(&tau1)
		merge2x1(&AL, &A11, &A21)
		// do an unblocked update to current block
		unblockedBuildQ(&AL, &tau1, &Wrk, 0)

		// zero upper part
		A01.SetIndexes(0.0)
		// --------------------------------------------------------
		continue3x3to2x2(
			&ATL, &ATR,
			&ABL, &ABR, &A00, &A11, &A22, A, pTOPLEFT)
		continue3x3to2x2(
			&TTL, &TTR,
			&TBL, &TBR, &T00, &T11, &T22, T, pTOPLEFT)
	}
	return err
}

/*
 * Generate an M-by-N real matrix Q with ortonormal columns, which is
 * defined as the product of k elementary reflectors and block reflector T
 *
 *    Q = H(1) H(2) . . . H(k)
 *
 * as returned by DecomposeQRT().
 *
 * Arguments:
 *  A     On entry, QR factorization as returned by DecomposeQRT() where the lower
 *        trapezoidal  part holds the elementary reflectors. On exit, the M-by-N
 *        matrix Q.
 *
 *  tau   The scalar factors of elementary reflectors as returned by   DecomposeQR()
 *
 *  W     Workspace, size A.Cols()-by-nb.
 *
 *  nb    Blocksize for blocked invocations. If nb == 0 unblocked algorith is used
 *
 * Compatible with lapack.DORGQR
 */
func BuildQ(A, tau, W *matrix.FloatMatrix, nb int) (*matrix.FloatMatrix, error) {
	var err error = nil
	if nb != 0 && W == nil {
		return nil, errors.New("workspace not defined")
	}
	// default is from LEFT
	if nb != 0 && (W.Cols() < nb || W.Rows() < A.Cols()) {
		return nil, errors.New("workspace too small")
	}

	if nb == 0 {
		w := matrix.FloatZeros(1, A.Cols())
		err = unblockedBuildQ(A, tau, w, 0)
	} else {
		var Wrk matrix.FloatMatrix
		Wrk.SubMatrixOf(W, 0, 0, A.Cols(), nb)
		err = blockedBuildQ(A, tau, &Wrk, nb)
	}
	return A, err
}

/*
 * Generate an M-by-N real matrix Q with ortonormal columns, which is
 * defined as the product of k elementary reflectors and block reflector T
 *
 *    Q = H(1) H(2) . . . H(k)
 *
 * generated using the compact WY representaion as returned by DecomposeQRT().
 *
 * Arguments:
 *  A     On entry, QR factorization as returned by DecomposeQRT() where the lower
 *        trapezoidal  part holds the elementary reflectors. On exit, the M-by-N
 *        matrix Q.
 *
 *  T     The block reflector computed from elementary reflectors as returned by
 *        DecomposeQRT() or computed from elementary reflectors and scalar coefficients
 *        by BuildT()
 *
 *  W     Workspace, size A.Cols()-by-nb.
 *
 *  nb    Blocksize for blocked invocations. If nb == 0 default value T.Cols()
 *        is used.
 *
 * Compatible with lapack.DORGQR
 */
func BuildQT(A, T, W *matrix.FloatMatrix, nb int) (*matrix.FloatMatrix, error) {
	var err error = nil
	if nb == 0 {
		nb = A.Cols()
	}
	// default is from LEFT
	if nb != 0 && (W.Cols() < nb || W.Rows() < A.Cols()) {
		return nil, errors.New("workspace too small")
	}

	var Wrk matrix.FloatMatrix
	Wrk.SubMatrixOf(W, 0, 0, A.Cols(), nb)
	err = blockedBuildQT(A, T, &Wrk, nb)
	return A, err
}

// Local Variables:
// tab-width: 4
// indent-tabs-mode: nil
// End:
