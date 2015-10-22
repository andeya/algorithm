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

/*
 * like LAPACK/dlafrt.f
 *
 * Build block reflector T from HH reflector stored in TriLU(A) and coefficients
 * in tau.
 *
 * Q = I - Y*T*Y.T; Householder H = I - tau*v*v.T
 *
 * T = | T  z |   z = -tau*T*Y.T*v
 *     | 0  c |   c = tau
 *
 * Q = H(1)H(2)...H(k) building forward here.
 */
func unblkQRBlockReflector(T, A, tau *matrix.FloatMatrix) {
	var ATL, ATR, ABL, ABR matrix.FloatMatrix
	var A00, a10, a11, A20, a21, A22 matrix.FloatMatrix
	var TTL, TTR, TBL, TBR matrix.FloatMatrix
	var T00, t01, T02, t11, t12, T22 matrix.FloatMatrix
	var tT, tB matrix.FloatMatrix
	var t0, tau1, t2 matrix.FloatMatrix

	partition2x2(
		&ATL, &ATR,
		&ABL, &ABR, A, 0, 0, pTOPLEFT)
	partition2x2(
		&TTL, &TTR,
		&TBL, &TBR, T, 0, 0, pTOPLEFT)
	partition2x1(
		&tT,
		&tB, tau, 0, pTOP)

	for ABR.Rows() > 0 && ABR.Cols() > 0 {
		repartition2x2to3x3(&ATL,
			&A00, nil, nil,
			&a10, &a11, nil,
			&A20, &a21, &A22, A, 1, pBOTTOMRIGHT)
		repartition2x2to3x3(&TTL,
			&T00, &t01, &T02,
			nil, &t11, &t12,
			nil, nil, &T22, T, 1, pBOTTOMRIGHT)
		repartition2x1to3x1(&tT,
			&t0,
			&tau1,
			&t2, tau, 1, pBOTTOM)
		// --------------------------------------------------

		// t11 := tau
		tauval := tau1.GetAt(0, 0)
		if tauval != 0.0 {
			t11.SetAt(0, 0, tauval)

			// t01 := a10.T + &A20.T*a21
			a10.CopyTo(&t01)
			MVMult(&t01, &A20, &a21, -tauval, -tauval, TRANSA)
			// t01 := T00*t01
			MVMultTrm(&t01, &T00, UPPER)
			//t01.Scale(-tauval)
		}

		// --------------------------------------------------
		continue3x3to2x2(
			&ATL, &ATR,
			&ABL, &ABR, &A00, &a11, &A22, A, pBOTTOMRIGHT)
		continue3x3to2x2(
			&TTL, &TTR,
			&TBL, &TBR, &T00, &t11, &T22, T, pBOTTOMRIGHT)
		continue3x1to2x1(
			&tT,
			&tB, &t0, &tau1, tau, pBOTTOM)
	}
}

/*
 * Unblocked QR decomposition with block reflector T.
 */
func unblockedQRT(A, T *matrix.FloatMatrix) {
	var ATL, ATR, ABL, ABR matrix.FloatMatrix
	var A00, a10, a11, a12, A20, a21, A22 matrix.FloatMatrix
	var TTL, TTR, TBL, TBR matrix.FloatMatrix
	var T00, t01, T02, t11, t12, T22 matrix.FloatMatrix

	//As.SubMatrixOf(A, 0, 0, mlen, nb)
	partition2x2(
		&ATL, &ATR,
		&ABL, &ABR, A, 0, 0, pTOPLEFT)
	partition2x2(
		&TTL, &TTR,
		&TBL, &TBR, T, 0, 0, pTOPLEFT)

	for ABR.Rows() > 0 && ABR.Cols() > 0 {
		repartition2x2to3x3(&ATL,
			&A00, nil, nil,
			&a10, &a11, &a12,
			&A20, &a21, &A22, A, 1, pBOTTOMRIGHT)
		repartition2x2to3x3(&TTL,
			&T00, &t01, &T02,
			nil, &t11, &t12,
			nil, nil, &T22, T, 1, pBOTTOMRIGHT)

		// ------------------------------------------------------

		computeHouseholder(&a11, &a21, &t11, LEFT)

		// H*[a12 A22].T
		applyHouseholder(&t11, &a21, &a12, &A22, LEFT)

		// update T
		tauval := t11.GetAt(0, 0)
		if tauval != 0.0 {
			// t01 := -tauval*(a10.T + &A20.T*a21)
			a10.CopyTo(&t01)
			MVMult(&t01, &A20, &a21, -tauval, -tauval, TRANSA)
			// t01 := T00*t01
			MVMultTrm(&t01, &T00, UPPER)
		}

		// ------------------------------------------------------
		continue3x3to2x2(
			&ATL, &ATR,
			&ABL, &ABR, &A00, &a11, &A22, A, pBOTTOMRIGHT)
		continue3x3to2x2(
			&TTL, &TTR,
			&TBL, &TBR, &T00, &t11, &T22, T, pBOTTOMRIGHT)
	}
}

/*
 * Unblocked QR decomposition. As implemented
 * in lapack.xGEQR2 subroutine.
 */
func unblockedQR(A, Tvec *matrix.FloatMatrix) {
	var ATL, ATR, ABL, ABR matrix.FloatMatrix
	var A00, a01, A02, a10, a11, a12, A20, a21, A22 matrix.FloatMatrix
	var tT, tB matrix.FloatMatrix
	var t0, tau1, t2 matrix.FloatMatrix

	//As.SubMatrixOf(A, 0, 0, mlen, nb)
	partition2x2(
		&ATL, &ATR,
		&ABL, &ABR, A, 0, 0, pTOPLEFT)
	partition2x1(
		&tT,
		&tB, Tvec, 0, pTOP)

	for ABR.Rows() > 0 && ABR.Cols() > 0 {
		repartition2x2to3x3(&ATL,
			&A00, &a01, &A02,
			&a10, &a11, &a12,
			&A20, &a21, &A22, A, 1, pBOTTOMRIGHT)
		repartition2x1to3x1(&tT,
			&t0,
			&tau1,
			&t2, Tvec, 1, pBOTTOM)

		// ------------------------------------------------------
		computeHouseholder(&a11, &a21, &tau1, LEFT)
		applyHouseholder(&tau1, &a21, &a12, &A22, LEFT)

		// ------------------------------------------------------
		continue3x3to2x2(
			&ATL, &ATR,
			&ABL, &ABR, &A00, &a11, &A22, A, pBOTTOMRIGHT)
		continue3x1to2x1(
			&tT,
			&tB, &t0, &tau1, Tvec, pBOTTOM)
	}
}

/*
 * Blocked QR decomposition with compact WY transform. As implemented
 * in lapack.xGEQRF subroutine.
 */
func blockedQR(A, Tvec, Twork, W *matrix.FloatMatrix, nb int) {
	var ATL, ATR, ABL, ABR, AL, AR matrix.FloatMatrix
	var A00, A01, A02, A10, A11, A12, A20, A21, A22 matrix.FloatMatrix
	var TT, TB matrix.FloatMatrix
	var t0, tau, t2, Tdiag, WT, WB, W0, W1, W2 matrix.FloatMatrix
	//var Twork, W *matrix.FloatMatrix

	Tdiag.DiagOf(Twork)

	partition2x2(
		&ATL, &ATR,
		&ABL, &ABR, A, 0, 0, pTOPLEFT)
	partition2x1(
		&TT,
		&TB, Tvec, 0, pTOP)
	partition2x1(
		&WT,
		&WB, W, 0, pTOP)

	for ABR.Rows() > 0 && ABR.Cols() > 0 {
		repartition2x2to3x3(&ATL,
			&A00, &A01, &A02,
			&A10, &A11, &A12,
			&A20, &A21, &A22, A, nb, pBOTTOMRIGHT)
		repartition2x1to3x1(&TT,
			&t0,
			&tau,
			&t2, Tvec, nb, pBOTTOM)
		repartition2x1to3x1(&WT,
			&W0,
			&W1,
			&W2, W, nb, pBOTTOM)
		partition1x2(
			&AL, &AR, &ABR, nb, pLEFT)

		// current block size
		cb, rb := A11.Size()
		if rb < cb {
			cb = rb
		}

		// --------------------------------------------------------

		// decompose left side AL == /A11\
		//                           \A21/
		unblockedQRT(&AL, Twork)

		// copy scaling from T diagonal to tau-vector
		Tdiag.CopyTo(&tau)

		// update A'tail i.e. A12 and A22 with (I - Y*T*Y.T).T * A'tail
		// compute: C - Y*(C.T*Y*T).T
		updateWithQT(&A12, &A22, &A11, &A21, Twork, &W2, cb, true)

		// --------------------------------------------------------
		continue3x3to2x2(
			&ATL, &ATR,
			&ABL, &ABR, &A00, &A11, &A22, A, pBOTTOMRIGHT)
		continue3x1to2x1(
			&TT,
			&TB, &t0, &tau, Tvec, pBOTTOM)
		continue3x1to2x1(
			&WT,
			&WB, &W0, &W1, W, pBOTTOM)
	}
}

// compute:
//      Q.T*C = (I -Y*T*Y.T).T*C ==  C - Y*(C.T*Y*T).T
// or
//      Q*C   = (I -Y*T*Y.T)*C   ==  C - Y*(C.T*Y*T.T).T
//
//
// where  C = /C1\   Y = /Y1\
//            \C2/       \Y2/
//
// C1 is nb*K, C2 is P*K, Y1 is nb*nb trilu, Y2 is P*nb, T is nb*nb
// W = K*nb
func updateWithQT(C1, C2, Y1, Y2, T, W *matrix.FloatMatrix, nb int, transpose bool) {

	/*
	   if transpose && W.Rows() != C1.Cols() {
	       panic(fmt.Sprintf("W.Rows [%d] != C1.Cols [%d]", W.Rows(), C1.Cols()))
	   } else if W.Rows() != C1.Rows() {
	       panic(fmt.Sprintf("W.Rows [%d] != C1.Rows [%d]", W.Rows(), C1.Rows()))
	   }
	*/

	// W = C1.T
	ScalePlus(W, C1, 0.0, 1.0, TRANSB)
	// W = C1.T*Y1
	MultTrm(W, Y1, 1.0, LOWER|UNIT|RIGHT)
	// W = W + C2.T*Y2
	Mult(W, C2, Y2, 1.0, 1.0, TRANSA)

	// --- here: W == C.T*Y ---
	tflags := UPPER | RIGHT
	if !transpose {
		tflags |= TRANSA
	}
	// W = W*T or W*T.T
	MultTrm(W, T, 1.0, Flags(tflags))

	// --- here: W == C.T*Y*T or C.T*Y*T.T ---

	// C2 = C2 - Y2*W.T
	Mult(C2, Y2, W, -1.0, 1.0, TRANSB)
	//  W = Y1*W.T ==> W.T = W*Y1.T
	MultTrm(W, Y1, 1.0, LOWER|UNIT|TRANSA|RIGHT)

	// C1 = C1 - W.T
	ScalePlus(C1, W, 1.0, -1.0, TRANSB)

	// --- here: C = (I - Y*T*Y.T).T * C ---
}

// compute:
//      C*Q.T = C*(I -Y*T*Y.T).T ==  C - C*Y*T.T*Y.T
// or
//      C*Q   = (I -Y*T*Y.T)*C   ==  C - C*Y*T*Y.T
//
//
// where  C = ( C1 C2 )   Y = ( Y1 )
//                            ( Y2 )
//
// C1 is K*nb, C2 is K*P, Y1 is nb*nb trilu, Y2 is P*nb, T is nb*nb
// W = K*nb
func updateWithQTRight(C1, C2, Y1, Y2, T, W *matrix.FloatMatrix, nb int, transpose bool) {

	// -- compute: W = C*Y = C1*Y1 + C2*Y2

	// W = C1
	ScalePlus(W, C1, 0.0, 1.0, NOTRANS)
	// W = C1*Y1
	MultTrm(W, Y1, 1.0, LOWER|UNIT|RIGHT)
	// W = W + C2*Y2
	Mult(W, C2, Y2, 1.0, 1.0, NOTRANS)

	// --- here: W == C*Y ---

	tflags := UPPER | RIGHT
	if transpose {
		tflags |= TRANSA
	}
	// W = W*T or W*T.T
	MultTrm(W, T, 1.0, Flags(tflags))

	// --- here: W == C*Y*T or C*Y*T.T ---

	// C2 = C2 - W*Y2.T
	Mult(C2, W, Y2, -1.0, 1.0, TRANSB)
	// C1 = C1 - W*Y1.T
	//  W = W*Y1
	MultTrm(W, Y1, 1.0, LOWER|UNIT|RIGHT|TRANSA)

	// C1 = C1 - W
	ScalePlus(C1, W, 1.0, -1.0, NOTRANS)

	// --- here: C = (I - Y*T*Y.T).T * C ---
}

func blockedQRT(A, T, W *matrix.FloatMatrix, nb int) {
	var ATL, ATR, ABL, ABR, AL, AR matrix.FloatMatrix
	var A00, A01, A02, A10, A11, A12, A20, A21, A22 matrix.FloatMatrix
	var WT, WB, W0, W1, W2 matrix.FloatMatrix
	var TTL, TTR, TBL, TBR matrix.FloatMatrix
	var T00, T01, T02, T11, T12, T22 matrix.FloatMatrix

	partition2x2(
		&ATL, &ATR,
		&ABL, &ABR, A, 0, 0, pTOPLEFT)
	partition2x2(
		&TTL, &TTR,
		&TBL, &TBR, T, 0, 0, pTOPLEFT)
	partition2x1(
		&WT,
		&WB, W, 0, pTOP)

	for ABR.Rows() > 0 && ABR.Cols() > 0 {
		repartition2x2to3x3(&ATL,
			&A00, &A01, &A02,
			&A10, &A11, &A12,
			&A20, &A21, &A22, A, nb, pBOTTOMRIGHT)
		repartition2x2to3x3(&TTL,
			&T00, &T01, &T02,
			nil, &T11, &T12,
			nil, nil, &T22, T, nb, pBOTTOMRIGHT)
		repartition2x1to3x1(&WT,
			&W0,
			&W1,
			&W2, W, nb, pBOTTOM)
		partition1x2(
			&AL, &AR, &ABR, nb, pLEFT)

		// current block size
		cb, rb := A11.Size()
		if rb < cb {
			cb = rb
		}

		// --------------------------------------------------------

		// decompose left side AL == /A11\
		//                           \A21/
		unblockedQRT(&AL, &T11)

		// update A'tail i.e. A12 and A22 with (I - Y*T*Y.T).T * A'tail
		// compute: Q*T.C == C - Y*(C.T*Y*T).T
		updateWithQT(&A12, &A22, &A11, &A21, &T11, &W2, cb, true)

		// update T01: T01 = -T00*Y1.T*Y2*T11
		//  Y1 = /A10\   Y2 = /A11\
		//       \A20/        \A21/
		//
		updateQRTReflector(&T01, &A10, &A20, &A11, &A21, &T00, &T11)

		// --------------------------------------------------------
		continue3x3to2x2(
			&ATL, &ATR,
			&ABL, &ABR, &A00, &A11, &A22, A, pBOTTOMRIGHT)
		continue3x3to2x2(
			&TTL, &TTR,
			&TBL, &TBR, &T00, &T11, &T22, T, pBOTTOMRIGHT)
		continue3x1to2x1(
			&WT,
			&WB, &W0, &W1, W, pBOTTOM)
	}
}

// update T: T = -T1*Y1.T*Y2*T2
//  Y1 = /Y10\   Y2 = /Y11\
//       \Y20/        \Y21/
//
//  T = -T1 * [Y10.T*Y11 + Y20.T*Y21]*T2
//
//  T1 is K*K tridiagonal upper matrix
//  T2 is nb*nb tridiagonal upper matrix
//  T  is K*nb block matrix
//  Y10 is nb*K block matrix
//  Y20 is M-K-nb*K block matrix
//  Y11 is nb*nb tridiagonal lower unit diagonal matrix
//  Y21 is M-K-nb*nb block matrix
//
func updateQRTReflector(T, Y10, Y20, Y11, Y21, T1, T2 *matrix.FloatMatrix) {
	// T = Y10.T
	if Y10.Cols() == 0 {
		return
	}
	// T = Y10.T
	ScalePlus(T, Y10, 0.0, 1.0, TRANSB)
	// T = Y10.T*Y11
	MultTrm(T, Y11, 1.0, LOWER|UNIT|RIGHT)
	// T = T + Y20.T*Y21
	Mult(T, Y20, Y21, 1.0, 1.0, TRANSA)
	// -- here: T == Y1.T*Y2

	// T = -T1*T
	MultTrm(T, T1, -1.0, UPPER)
	// T = T*T2
	MultTrm(T, T2, 1.0, UPPER|RIGHT)
}

/*
 * Compute QR factorization of a M-by-N matrix A: A = Q * R.
 *
 * Arguments:
 *  A   On entry, the M-by-N matrix A. On exit, the elements on and above
 *      the diagonal contain the min(M,N)-by-N upper trapezoidal matrix R.
 *      The elements below the diagonal with the column vector 'tau', represent
 *      the ortogonal matrix Q as product of elementary reflectors.
 *
 * tau  On exit, the scalar factors of the elemenentary reflectors.
 *
 * W    Workspace, N-by-nb matrix used for work space in blocked invocations.
 *
 * nb   The block size used in blocked invocations. If nb is zero on N <= nb
 *      unblocked algorithm is used.
 *
 * Returns:
 *      Decomposed matrix A and error indicator.
 *
 * DecomposeQR is compatible with lapack.DGEQRF
 */
func DecomposeQR(A, tau, W *matrix.FloatMatrix, nb int) (*matrix.FloatMatrix, error) {
	var err error = nil

	if nb == 0 || A.Cols() <= nb {
		unblockedQR(A, tau)
	} else {
		Twork := matrix.FloatZeros(nb, nb)
		if W == nil {
			W = matrix.FloatZeros(A.Cols(), nb)
		} else if W.Cols() < nb || W.Rows() < A.Cols() {
			return nil, errors.New("work space too small")
		}
		var Wrk matrix.FloatMatrix
		Wrk.SubMatrixOf(W, 0, 0, A.Cols(), nb)
		blockedQR(A, tau, Twork, &Wrk, nb)
	}
	return A, err
}

/*
 * Compute QR factorization of a M-by-N matrix A using compact WY transformation: A = Q * R,
 * where Q = I - Y*T*Y.T, T is block reflector and Y holds elementary reflectors as lower
 * trapezoidal matrix saved below diagonal elements of the matrix A.
 *
 * Arguments:
 *  A   On entry, the M-by-N matrix A. On exit, the elements on and above
 *      the diagonal contain the min(M,N)-by-N upper trapezoidal matrix R.
 *      The elements below the diagonal with the matrix 'T', represent
 *      the ortogonal matrix Q as product of elementary reflectors.
 *
 * T    On exit, the block reflector which, together with trilu(A) represent
 *      the ortogonal matrix Q as Q = I - Y*T*Y.T where Y = trilu(A).
 *
 * W    Workspace, N-by-nb matrix used for work space in blocked invocations.
 *
 * nb   The block size used in blocked invocations. If nb is zero on N <= nb
 *      unblocked algorithm is used.
 *
 * Returns:
 *      Decomposed matrix A and error indicator.
 *
 * DecomposeQRT is compatible with lapack.DGEQRT
 */
func DecomposeQRT(A, T, W *matrix.FloatMatrix, nb int) (*matrix.FloatMatrix, error) {
	var err error = nil
	if nb == 0 || A.Cols() <= nb {
		unblockedQRT(A, T)
	} else {
		if W == nil {
			W = matrix.FloatZeros(A.Cols(), nb)
		} else if W.Cols() < nb || W.Rows() < A.Cols() {
			return nil, errors.New("work space too small")
		}
		var Wrk matrix.FloatMatrix
		Wrk.SubMatrixOf(W, 0, 0, A.Cols(), nb)
		blockedQRT(A, T, &Wrk, nb)
	}
	return A, err
}

/*
 * Build block reflector T from HH elementary reflectors stored in TriLU(A) and
 * scalar factors in tau.
 *
 * Q = I - Y*T*Y.T; Householder H = I - tau*v*v.T
 *
 * T = | T  z |   z = -tau*T*Y.T*v
 *     | 0  c |   c = tau
 *
 * Compatible with lapack.DLAFRT
 */
func BuildT(T, A, tau *matrix.FloatMatrix) (*matrix.FloatMatrix, error) {
	var err error = nil

	if T.Cols() < A.Cols() || T.Rows() < A.Cols() {
		return nil, errors.New("reflector matrix T too small")
	}

	unblkQRBlockReflector(T, A, tau)
	return T, err
}

// Local Variables:
// tab-width: 4
// indent-tabs-mode: nil
// End:
