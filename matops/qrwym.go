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
 * Unblocked algorith for computing C = Q.T*C and C = Q*C.
 *
 * Q = H(1)H(2)...H(k) where elementary reflectors H(i) are stored on i'th column
 * below diagonal in A.
 *
 * Progressing A from top-left to bottom-right i.e from smaller column numbers
 * to larger, produces H(k)...H(2)H(1) == Q.T. and C = Q.T*C
 *
 * Progressing from bottom-right to top-left produces H(1)H(2)...H(k) == Q and C = Q*C
 */
func unblockedMultQLeft(C, A, tau, w *matrix.FloatMatrix, flags Flags) {
	var ATL, ATR, ABL, ABR matrix.FloatMatrix
	var A00, a10, a11, A20, a21, A22 matrix.FloatMatrix
	var CT, CB, C0, c1t, C2 matrix.FloatMatrix
	var tT, tB matrix.FloatMatrix
	var t0, tau1, t2 matrix.FloatMatrix

	var Aref *matrix.FloatMatrix
	var pAdir, pAstart, pDir, pStart pDirection
	var mb int

	// partitioning start and direction
	if flags&TRANS != 0 {
		// from top-left to bottom-right to produce transposed sequence (Q.T*C)
		pAstart = pTOPLEFT
		pAdir = pBOTTOMRIGHT
		pStart = pTOP
		pDir = pBOTTOM
		mb = 0
		Aref = &ABR
	} else {
		// from bottom-right to top-left to produce normal sequence (Q*C)
		pAstart = pBOTTOMRIGHT
		pAdir = pTOPLEFT
		pStart = pBOTTOM
		pDir = pTOP
		mb = A.Rows() - A.Cols()
		Aref = &ATL
	}

	partition2x2(
		&ATL, &ATR,
		&ABL, &ABR, A, mb, 0, pAstart)
	partition2x1(
		&CT,
		&CB, C, mb, pStart)
	partition2x1(
		&tT,
		&tB, tau, 0, pStart)

	for Aref.Rows() > 0 && Aref.Cols() > 0 {
		repartition2x2to3x3(&ATL,
			&A00, nil, nil,
			&a10, &a11, nil,
			&A20, &a21, &A22, A, 1, pAdir)
		repartition2x1to3x1(&CT,
			&C0,
			&c1t,
			&C2, C, 1, pDir)
		repartition2x1to3x1(&tT,
			&t0,
			&tau1,
			&t2, tau, 1, pDir)

		// --------------------------------------------------------

		applyHHTo2x1(&tau1, &a21, &c1t, &C2, w, LEFT)

		// --------------------------------------------------------
		continue3x3to2x2(
			&ATL, &ATR,
			&ABL, &ABR, &A00, &a11, &A22, A, pAdir)
		continue3x1to2x1(
			&CT,
			&CB, &C0, &c1t, C, pDir)
		continue3x1to2x1(
			&tT,
			&tB, &t0, &tau1, tau, pDir)
	}
}

/*
 * Unblocked algorith for computing C = C*Q.T and C = C*Q.
 *
 * Q = H(1)H(2)...H(k) where elementary reflectors H(i) are stored on i'th column
 * below diagonal in A.
 *
 *     Q.T = (H1(1)*H(2)*...*H(k)).T
 *         = H(k).T*...*H(2).T*H(1).T
 *         = H(k)...H(2)H(1)
 *
 * Progressing A from top-left to bottom-right i.e from smaller column numbers
 * to larger, produces C*H(1)H(2)...H(k) == C*Q.
 *
 * Progressing from bottom-right to top-left produces C*H(k)...H(2)H(1) == C*Q.T.
 */
func unblockedMultQRight(C, A, tau, w *matrix.FloatMatrix, flags Flags) {
	var ATL, ATR, ABL, ABR matrix.FloatMatrix
	var A00, a10, a11, A20, a21, A22 matrix.FloatMatrix
	var CL, CR, C0, c1, C2 matrix.FloatMatrix
	var tT, tB matrix.FloatMatrix
	var t0, tau1, t2 matrix.FloatMatrix

	var Aref *matrix.FloatMatrix
	var pAdir, pAstart, pDir, pStart, pCstart, pCdir pDirection
	var cb, mb int

	// partitioning start and direction
	if flags&TRANS != 0 {
		// from bottom-right to top-left to produce transpose sequence (C*Q.T)
		pAstart = pBOTTOMRIGHT
		pAdir = pTOPLEFT
		pStart = pBOTTOM
		pDir = pTOP
		pCstart = pRIGHT
		pCdir = pLEFT
		mb = A.Rows() - A.Cols()
		cb = C.Cols() - A.Cols()
		Aref = &ATL
	} else {
		// from top-left to bottom-right to produce normal sequence (C*Q)
		pAstart = pTOPLEFT
		pAdir = pBOTTOMRIGHT
		pStart = pTOP
		pDir = pBOTTOM
		pCstart = pLEFT
		pCdir = pRIGHT
		mb = 0
		cb = 0
		Aref = &ABR
	}

	partition2x2(
		&ATL, &ATR,
		&ABL, &ABR, A, mb, 0, pAstart)
	partition1x2(
		&CL, &CR, C, cb, pCstart)
	partition2x1(
		&tT,
		&tB, tau, 0, pStart)

	for Aref.Rows() > 0 && Aref.Cols() > 0 {
		repartition2x2to3x3(&ATL,
			&A00, nil, nil,
			&a10, &a11, nil,
			&A20, &a21, &A22, A, 1, pAdir)
		repartition1x2to1x3(&CL,
			&C0, &c1, &C2, C, 1, pCdir)
		repartition2x1to3x1(&tT,
			&t0,
			&tau1,
			&t2, tau, 1, pDir)

		// --------------------------------------------------------

		applyHHTo2x1(&tau1, &a21, &c1, &C2, w, RIGHT)

		// --------------------------------------------------------
		continue3x3to2x2(
			&ATL, &ATR,
			&ABL, &ABR, &A00, &a11, &A22, A, pAdir)
		continue1x3to1x2(
			&CL, &CR, &C0, &c1, C, pCdir)
		continue3x1to2x1(
			&tT,
			&tB, &t0, &tau1, tau, pDir)
	}
}

/*
 * Blocked version for computing C = Q*C and C = Q.T*C from elementary reflectors
 * and scalar coefficients.
 *
 * Elementary reflectors and scalar coefficients are used to build block reflector T.
 * Matrix C is updated by applying block reflector T using compact WY algorithm.
 */
func blockedMultQLeft(C, A, tau, W *matrix.FloatMatrix, nb int, flags Flags) {
	var ATL, ATR, ABL, ABR, AL matrix.FloatMatrix
	var A00, A10, A11, A20, A21, A22 matrix.FloatMatrix
	var CT, CB, C0, C1, C2 matrix.FloatMatrix
	var tT, tB matrix.FloatMatrix
	var t0, tau1, t2 matrix.FloatMatrix
	var Wrk, Tw matrix.FloatMatrix

	var Aref *matrix.FloatMatrix
	var pAdir, pAstart, pDir, pStart pDirection
	var bsz, mb int

	// partitioning start and direction
	if flags&TRANS != 0 || nb == A.Cols() {
		// from top-left to bottom-right to produce transposed sequence (Q.T*C)
		pAstart = pTOPLEFT
		pAdir = pBOTTOMRIGHT
		pStart = pTOP
		pDir = pBOTTOM
		mb = 0
		Aref = &ABR
	} else {
		// from bottom-right to top-left to produce normal sequence (Q*C)
		pAstart = pBOTTOMRIGHT
		pAdir = pTOPLEFT
		pStart = pBOTTOM
		pDir = pTOP
		mb = A.Rows() - A.Cols()
		Aref = &ATL
	}

	Twork := matrix.FloatZeros(nb, nb)

	partition2x2(
		&ATL, &ATR,
		&ABL, &ABR, A, mb, 0, pAstart)
	partition2x1(
		&CT,
		&CB, C, mb, pStart)
	partition2x1(
		&tT,
		&tB, tau, 0, pStart)

	transpose := flags&TRANS != 0

	for Aref.Rows() > 0 && Aref.Cols() > 0 {
		repartition2x2to3x3(&ATL,
			&A00, nil, nil,
			&A10, &A11, nil,
			&A20, &A21, &A22, A, nb, pAdir)
		repartition2x1to3x1(&tT,
			&t0,
			&tau1,
			&t2, tau, nb, pDir)
		bsz = A11.Cols()
		repartition2x1to3x1(&CT,
			&C0,
			&C1,
			&C2, C, bsz, pDir)

		// --------------------------------------------------------

		// build block reflector from current block
		merge2x1(&AL, &A11, &A21)
		Tw.SubMatrixOf(Twork, 0, 0, bsz, bsz)
		unblkQRBlockReflector(&Tw, &AL, &tau1)

		// compute: Q*T.C == C - Y*(C.T*Y*T).T  transpose == true
		//          Q*C   == C - C*Y*T*Y.T      transpose == false
		Wrk.SubMatrixOf(W, 0, 0, C1.Cols(), bsz)
		updateWithQT(&C1, &C2, &A11, &A21, &Tw, &Wrk, nb, transpose)

		// --------------------------------------------------------
		continue3x3to2x2(
			&ATL, &ATR,
			&ABL, &ABR, &A00, &A11, &A22, A, pAdir)
		continue3x1to2x1(
			&CT,
			&CB, &C0, &C1, C, pDir)
		continue3x1to2x1(
			&tT,
			&tB, &t0, &tau1, tau, pDir)
	}

}

/*
 * Blocked version for computing C = C*Q and C = C*Q.T from elementary reflectors
 * and scalar coefficients.
 *
 * Elementary reflectors and scalar coefficients are used to build block reflector T.
 * Matrix C is updated by applying block reflector T using compact WY algorithm.
 */
func blockedMultQRight(C, A, tau, W *matrix.FloatMatrix, nb int, flags Flags) {
	var ATL, ATR, ABL, ABR, AL matrix.FloatMatrix
	var A00, A10, A11, A20, A21, A22 matrix.FloatMatrix
	var CL, CR, C0, C1, C2 matrix.FloatMatrix
	var tT, tB matrix.FloatMatrix
	var t0, tau1, t2 matrix.FloatMatrix
	var Wrk, Tw matrix.FloatMatrix

	var Aref *matrix.FloatMatrix
	var pAdir, pAstart, pDir, pStart, pCstart, pCdir pDirection
	var bsz, cb, mb int

	// partitioning start and direction
	if flags&TRANS != 0 {
		// from bottom-right to top-left to produce transpose sequence (C*Q.T)
		pAstart = pBOTTOMRIGHT
		pAdir = pTOPLEFT
		pStart = pBOTTOM
		pDir = pTOP
		pCstart = pRIGHT
		pCdir = pLEFT
		mb = A.Rows() - A.Cols()
		cb = C.Cols() - A.Cols()
		Aref = &ATL
	} else {
		// from top-left to bottom-right to produce normal sequence (C*Q)
		pAstart = pTOPLEFT
		pAdir = pBOTTOMRIGHT
		pStart = pTOP
		pDir = pBOTTOM
		pCstart = pLEFT
		pCdir = pRIGHT
		mb = 0
		cb = 0
		Aref = &ABR
	}

	Twork := matrix.FloatZeros(nb, nb)

	partition2x2(
		&ATL, &ATR,
		&ABL, &ABR, A, mb, 0, pAstart)
	partition1x2(
		&CL, &CR, C, cb, pCstart)
	partition2x1(
		&tT,
		&tB, tau, 0, pStart)

	transpose := flags&TRANS != 0

	for Aref.Rows() > 0 && Aref.Cols() > 0 {
		repartition2x2to3x3(&ATL,
			&A00, nil, nil,
			&A10, &A11, nil,
			&A20, &A21, &A22, A, nb, pAdir)
		repartition2x1to3x1(&tT,
			&t0,
			&tau1,
			&t2, tau, nb, pDir)

		bsz = A11.Cols() // C1 block size must match A11
		repartition1x2to1x3(&CL,
			&C0, &C1, &C2, C, bsz, pCdir)

		// --------------------------------------------------------

		// build block reflector from current block
		merge2x1(&AL, &A11, &A21)
		Tw.SubMatrixOf(Twork, 0, 0, bsz, bsz)
		unblkQRBlockReflector(&Tw, &AL, &tau1)

		// compute: C*Q.T == C - C*(Y*T*Y.T).T = C - C*Y*T.T*Y.T
		//          C*Q   == C - C*Y*T*Y.T
		Wrk.SubMatrixOf(W, 0, 0, C1.Rows(), bsz)
		updateWithQTRight(&C1, &C2, &A11, &A21, &Tw, &Wrk, nb, transpose)

		// --------------------------------------------------------
		continue3x3to2x2(
			&ATL, &ATR,
			&ABL, &ABR, &A00, &A11, &A22, A, pAdir)
		continue1x3to1x2(
			&CL, &CR, &C0, &C1, C, pCdir)
		continue3x1to2x1(
			&tT,
			&tB, &t0, &tau1, tau, pDir)
	}

}

/*
 * Blocked version for computing C = Q*C and C = Q.T*C with block reflector.
 *
 */
func blockedMultQTLeft(C, A, T, W *matrix.FloatMatrix, nb int, flags Flags) {
	var ATL, ATR, ABL, ABR matrix.FloatMatrix
	var A00, A10, A11, A20, A21, A22 matrix.FloatMatrix
	var CT, CB, C0, C1, C2 matrix.FloatMatrix
	var TTL, TTR, TBL, TBR matrix.FloatMatrix
	var T00, T01, T02, T11, T12, T22 matrix.FloatMatrix

	var Aref *matrix.FloatMatrix
	var pAdir, pAstart, pCdir, pCstart pDirection
	var bsz, mb int

	// partitioning start and direction
	if flags&TRANS != 0 {
		// from top-left to bottom-right to produce transposed sequence (Q.T*C)
		pAstart = pTOPLEFT
		pAdir = pBOTTOMRIGHT
		pCstart = pTOP
		pCdir = pBOTTOM
		mb = 0
		Aref = &ABR
	} else {
		// from bottom-right to top-left to produce normal sequence (Q*C)
		pAstart = pBOTTOMRIGHT
		pAdir = pTOPLEFT
		pCstart = pBOTTOM
		pCdir = pTOP
		mb = A.Rows() - A.Cols()
		Aref = &ATL
	}

	partition2x2(
		&ATL, &ATR,
		&ABL, &ABR, A, mb, 0, pAstart)
	partition2x2(
		&TTL, &TTR,
		&TBL, &TBR, T, 0, 0, pAstart)
	partition2x1(
		&CT,
		&CB, C, mb, pCstart)

	transpose := flags&TRANS != 0

	for Aref.Rows() > 0 && Aref.Cols() > 0 {
		repartition2x2to3x3(&ATL,
			&A00, nil, nil,
			&A10, &A11, nil,
			&A20, &A21, &A22, A, nb, pAdir)
		repartition2x2to3x3(&TTL,
			&T00, &T01, &T02,
			nil, &T11, &T12,
			nil, nil, &T22, T, nb, pAdir)

		bsz = A11.Cols() // must match A11 block size
		repartition2x1to3x1(&CT,
			&C0,
			&C1,
			&C2, C, bsz, pCdir)

		// --------------------------------------------------------

		// compute: Q.T*C == C - Y*(C.T*Y*T).T  transpose == true
		//          Q*C   == C - C*Y*T*Y.T      transpose == false

		var Wrk matrix.FloatMatrix
		Wrk.SubMatrixOf(W, 0, 0, C1.Cols(), bsz)
		updateWithQT(&C1, &C2, &A11, &A21, &T11, &Wrk, nb, transpose)

		// --------------------------------------------------------
		continue3x3to2x2(
			&ATL, &ATR,
			&ABL, &ABR, &A00, &A11, &A22, A, pAdir)
		continue3x3to2x2(
			&TTL, &TTR,
			&TBL, &TBR, &T00, &T11, &T22, T, pAdir)
		continue3x1to2x1(
			&CT,
			&CB, &C0, &C1, C, pCdir)
	}

}

/*
 * Blocked version for computing C = C*Q and C = C*Q.T with block reflector.
 *
 */
func blockedMultQTRight(C, A, T, W *matrix.FloatMatrix, nb int, flags Flags) {
	var ATL, ATR, ABL, ABR matrix.FloatMatrix
	var A00, A10, A11, A20, A21, A22 matrix.FloatMatrix
	var CL, CR, C0, C1, C2 matrix.FloatMatrix
	var TTL, TTR, TBL, TBR matrix.FloatMatrix
	var T00, T01, T02, T11, T12, T22 matrix.FloatMatrix

	var Aref *matrix.FloatMatrix
	var pAdir, pAstart, pCstart, pCdir pDirection
	var bsz, cb, mb int

	// partitioning start and direction
	if flags&TRANS != 0 {
		// from bottom-right to top-left to produce transpose sequence (C*Q.T)
		pAstart = pBOTTOMRIGHT
		pAdir = pTOPLEFT
		pCstart = pRIGHT
		pCdir = pLEFT
		mb = A.Rows() - A.Cols()
		cb = C.Cols() - A.Cols()
		Aref = &ATL
	} else {
		// from top-left to bottom-right to produce normal sequence (C*Q)
		pAstart = pTOPLEFT
		pAdir = pBOTTOMRIGHT
		pCstart = pLEFT
		pCdir = pRIGHT
		mb = 0
		cb = 0
		Aref = &ABR
	}

	partition2x2(
		&ATL, &ATR,
		&ABL, &ABR, A, mb, 0, pAstart)
	partition2x2(
		&TTL, &TTR,
		&TBL, &TBR, T, 0, 0, pAstart)
	partition1x2(
		&CL, &CR, C, cb, pCstart)

	transpose := flags&TRANS != 0

	for Aref.Rows() > 0 && Aref.Cols() > 0 {
		repartition2x2to3x3(&ATL,
			&A00, nil, nil,
			&A10, &A11, nil,
			&A20, &A21, &A22, A, nb, pAdir)
		repartition2x2to3x3(&TTL,
			&T00, &T01, &T02,
			nil, &T11, &T12,
			nil, nil, &T22, T, nb, pAdir)

		bsz = A11.Cols()
		repartition1x2to1x3(&CL,
			&C0, &C1, &C2, C, bsz, pCdir)

		// --------------------------------------------------------

		// compute: C*Q.T == C - C*Y*T.T*Y.T   transpose == true
		//          C*Q   == C - C*Y*T*Y.T     transpose == false

		var Wrk matrix.FloatMatrix
		Wrk.SubMatrixOf(W, 0, 0, C1.Rows(), bsz)
		updateWithQTRight(&C1, &C2, &A11, &A21, &T11, &Wrk, nb, transpose)

		// --------------------------------------------------------
		continue3x3to2x2(
			&ATL, &ATR,
			&ABL, &ABR, &A00, &A11, &A22, A, pAdir)
		continue3x3to2x2(
			&TTL, &TTR,
			&TBL, &TBR, &T00, &T11, &T22, T, pAdir)
		continue1x3to1x2(
			&CL, &CR, &C0, &C1, C, pCdir)
	}

}

/*
 * Multiply and replace C with Q*C or Q.T*C where Q is a real orthogonal matrix
 * defined as the product of k elementary reflectors.
 *
 *    Q = H(1) H(2) . . . H(k)
 *
 * as returned by DecomposeQR().
 *
 * Arguments:
 *  C     On entry, the M-by-N matrix C. On exit C is overwritten by Q*C or Q.T*C.
 *
 *  A     QR factorization as returne by DecomposeQR() where the lower trapezoidal
 *        part holds the elementary reflectors.
 *
 *  tau   The scalar factors of the elementary reflectors.
 *
 *  W     Workspace, used for blocked invocations. Size C.Cols()-by-nb.
 *
 *  nb    Blocksize for blocked invocations. If C.Cols() <= nb unblocked algorithm
 *        is used.
 *
 *  flags Indicators. Valid indicators LEFT, RIGHT, TRANS, NOTRANS
 *
 * Compatible with lapack.DORMQR
 */
func MultQ(C, A, tau, W *matrix.FloatMatrix, flags Flags, nb int) error {
	var err error = nil
	if nb != 0 && W == nil {
		return errors.New("workspace not defined")
	}
	if flags&RIGHT != 0 {
		// from right; C*A or C*A.T
		if C.Cols() != A.Rows() {
			return errors.New("C*Q: C.Cols != A.Rows")
		}
		if nb != 0 && (W.Cols() < nb || W.Rows() < C.Rows()) {
			return errors.New("workspace too small")
		}
	} else {
		// default is from LEFT; A*C or A.T*C
		/*
		   if C.Rows() != A.Rows() {
		       return errors.New("Q*C: C.Rows != A.Rows")
		   }
		*/
		if nb != 0 && (W.Cols() < nb || W.Rows() < C.Cols()) {
			return errors.New("workspace too small")
		}
	}
	if nb == 0 {
		if flags&RIGHT != 0 {
			w := matrix.FloatZeros(C.Rows(), 1)
			unblockedMultQRight(C, A, tau, w, flags)
		} else {
			w := matrix.FloatZeros(1, C.Cols())
			unblockedMultQLeft(C, A, tau, w, flags)
		}
	} else {
		var Wrk matrix.FloatMatrix
		if flags&RIGHT != 0 {
			Wrk.SubMatrixOf(W, 0, 0, C.Rows(), nb)
			blockedMultQRight(C, A, tau, &Wrk, nb, flags)
		} else {
			Wrk.SubMatrixOf(W, 0, 0, C.Cols(), nb)
			blockedMultQLeft(C, A, tau, &Wrk, nb, flags)
		}
	}
	return err
}

/*
 * Multiply and replace C with Q*C or Q.T*C where Q is a real orthogonal matrix
 * defined as the product of k elementary reflectors and block reflector T
 *
 *    Q = H(1) H(2) . . . H(k)
 *
 * as returned by DecomposeQRT().
 *
 * Arguments:
 *  C     On entry, the M-by-N matrix C. On exit C is overwritten by Q*C or Q.T*C.
 *
 *  A     QR factorization as returned by DecomposeQRT() where the lower trapezoidal
 *        part holds the elementary reflectors.
 *
 *  T     The block reflector computed from elementary reflectors as returned by
 *        DecomposeQRT() or computed from elementary reflectors and scalar coefficients
 *        by BuildT()
 *
 *  W     Workspace, size C.Cols()-by-nb or C.Rows()-by-nb
 *
 *  nb    Blocksize for blocked invocations. If nb == 0 default value T.Cols()
 *        is used.
 *
 *  flags Indicators. Valid indicators LEFT, RIGHT, TRANS, NOTRANS
 *
 * Compatible with lapack.DGEMQRT
 */
func MultQT(C, A, T, W *matrix.FloatMatrix, flags Flags, nb int) error {
	var err error = nil
	if nb == 0 {
		nb = T.Cols()
	}
	if W == nil {
		return errors.New("workspace not defined")
	}
	if flags&RIGHT != 0 {
		// from right; C*A or C*A.T
		if C.Cols() != A.Rows() {
			return errors.New("C*Q: C.Cols != A.Rows")
		}
		if W.Cols() < nb || W.Rows() < C.Rows() {
			return errors.New("workspace too small")
		}
	} else {
		// default is from LEFT; A*C or A.T*C
		/*
		   if C.Rows() != A.Rows() {
		       return errors.New("Q*C: C.Rows != A.Rows")
		   }
		*/
		if W.Cols() < nb || W.Rows() < C.Cols() {
			return errors.New("workspace too small")
		}
	}

	var Wrk matrix.FloatMatrix
	if flags&RIGHT != 0 {
		Wrk.SubMatrixOf(W, 0, 0, C.Rows(), nb)
		blockedMultQTRight(C, A, T, &Wrk, nb, flags)
	} else {
		Wrk.SubMatrixOf(W, 0, 0, C.Cols(), nb)
		blockedMultQTLeft(C, A, T, &Wrk, nb, flags)
	}
	return err
}

/*
 * Solve a system of linear equations A*X = B with general M-by-N
 * matrix A using the QR factorization computed by DecomposeQR().
 *
 * If flags&TRANS != 0:
 *   find the minimum norm solution of an overdetermined system A.T * X = B.
 *   i.e min ||X|| s.t A.T*X = B
 *
 * Otherwise:
 *   find the least squares solution of an overdetermined system, i.e.,
 *   solve the least squares problem: min || B - A*X ||.
 *
 * Arguments:
 *  B    On entry, the right hand side N-by-P matrix B. On exit, the solution matrix X.
 *
 *  A    The elements on and above the diagonal contain the min(M,N)-by-N upper
 *       trapezoidal matrix R. The elements below the diagonal with the vector 'tau',
 *       represent the ortogonal matrix Q as product of elementary reflectors.
 *       Matrix A and T are as returned by DecomposeQR()
 *
 *  tau  The vector of N scalar coefficients that together with trilu(A) define
 *       the ortogonal matrix Q as Q = H(1)H(2)...H(N)
 *
 *  W    Workspace, P-by-nb matrix used for work space in blocked invocations.
 *
 *  flags Indicator flag
 *
 *  nb   The block size used in blocked invocations. If nb is zero or P < nb
 *       unblocked algorithm is used.
 *
 * Compatible with lapack.GELS (the m >= n part)
 */
func SolveQR(B, A, tau, W *matrix.FloatMatrix, flags Flags, nb int) error {
	var err error = nil
	var R, BT matrix.FloatMatrix
	if flags&TRANS != 0 {
		// Solve overdetermined system A.T*X = B

		// B' = R.-1*B
		A.SubMatrix(&R, 0, 0, A.Cols(), A.Cols())
		B.SubMatrix(&BT, 0, 0, A.Cols(), B.Cols())
		err = SolveTrm(&BT, &R, 1.0, LEFT|UPPER|TRANSA)

		// Clear bottom part of B
		B.SubMatrixOf(&BT, A.Cols(), 0)
		BT.SetIndexes(0.0)

		// X = Q*B'
		err = MultQ(B, A, tau, W, LEFT, nb)
	} else {
		// solve least square problem min ||A*X - B||

		// B' = Q.T*B
		err = MultQ(B, A, tau, W, LEFT|TRANS, nb)
		if err != nil {
			return err
		}

		// X = R.-1*B'
		A.SubMatrix(&R, 0, 0, A.Cols(), A.Cols())
		B.SubMatrix(&BT, 0, 0, A.Cols(), B.Cols())
		err = SolveTrm(&BT, &R, 1.0, LEFT|UPPER)

	}
	return err
}

/*
 * Solve a system of linear equations A*X = B with general M-by-N
 * matrix A using the QR factorization computed by DecomposeQRT().
 *
 * If flags&TRANS != 0:
 *   find the minimum norm solution of an overdetermined system A.T * X = B.
 *   i.e min ||X|| s.t A.T*X = B
 *
 * Otherwise:
 *   find the least squares solution of an overdetermined system, i.e.,
 *   solve the least squares problem: min || B - A*X ||.
 *
 * Arguments:
 *  B     On entry, the right hand side N-by-P matrix B. On exit, the solution matrix X.
 *
 *  A     The elements on and above the diagonal contain the min(M,N)-by-N upper
 *        trapezoidal matrix R. The elements below the diagonal with the matrix 'T',
 *        represent the ortogonal matrix Q as product of elementary reflectors.
 *        Matrix A and T are as returned by DecomposeQRT()
 *
 *  T     The N-by-N block reflector which, together with trilu(A) represent
 *        the ortogonal matrix Q as Q = I - Y*T*Y.T where Y = trilu(A).
 *
 *  W     Workspace, P-by-nb matrix used for work space in blocked invocations.
 *
 *  flags Indicator flag
 *
 *  nb    The block size used in blocked invocations. If nb is zero default
 *        value N is used.
 *
 * Compatible with lapack.GELS (the m >= n part)
 */
func SolveQRT(B, A, T, W *matrix.FloatMatrix, flags Flags, nb int) error {
	var err error = nil
	var R, BT matrix.FloatMatrix
	if flags&TRANS != 0 {
		// Solve overdetermined system A.T*X = B

		// B' = R.-1*B
		A.SubMatrix(&R, 0, 0, A.Cols(), A.Cols())
		B.SubMatrix(&BT, 0, 0, A.Cols(), B.Cols())
		err = SolveTrm(&BT, &R, 1.0, LEFT|UPPER|TRANSA)

		// Clear bottom part of B
		B.SubMatrix(&BT, A.Cols(), 0)
		BT.SetIndexes(0.0)

		// X = Q*B'
		err = MultQT(B, A, T, W, LEFT, nb)
	} else {
		// solve least square problem min ||A*X - B||

		// B' = Q.T*B
		err = MultQT(B, A, T, W, LEFT|TRANS, nb)
		if err != nil {
			return err
		}

		// X = R.-1*B'
		A.SubMatrix(&R, 0, 0, A.Cols(), A.Cols())
		B.SubMatrix(&BT, 0, 0, A.Cols(), B.Cols())
		err = SolveTrm(&BT, &R, 1.0, LEFT|UPPER)
	}
	return err
}

// Local Variables:
// tab-width: 4
// indent-tabs-mode: nil
// End:
