// Copyright (c) Harri Rautila, 2012,2013

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
 *  ( a11  a12 )   ( 1   0   )( d1  0   )( l  l21.t )
 *  ( a21  A22 )   ( l21 L22 )(  0  A22 )( 0  L22.t )
 *
 *   a11  =   d1
 *   a21  =   l21*d1                       => l21 = a21/d1
 *   A22  =   l21*d1*l21.t + L22*D2*L22.t  => L22 = A22 - l21*d1*l21t
 */
func unblkLowerLDL(A *matrix.FloatMatrix, p *pPivots) (err error) {
	var ATL, ATR, ABL, ABR matrix.FloatMatrix
	var A00, a10, a11, A20, a21, A22, acol matrix.FloatMatrix
	var pT, pB, p0, p1, p2 pPivots

	err = nil
	partition2x2(
		&ATL, &ATR,
		&ABL, &ABR, A, 0, 0, pTOPLEFT)
	partitionPivot2x1(
		&pT,
		&pB, p, 0, pTOP)

	for ATL.Rows() < A.Rows() {
		repartition2x2to3x3(&ATL,
			&A00, nil, nil,
			&a10, &a11, nil,
			&A20, &a21, &A22, A, 1, pBOTTOMRIGHT)
		repartPivot2x1to3x1(&pT,
			&p0, &p1, &p2 /**/, p, 1, pBOTTOM)

		// --------------------------------------------------------

		ABR.Diag(&acol)
		//merge2x1(&acol, &a11, &a21)
		imax := IAMax(&acol)
		//fmt.Printf("imax=%d, val=%e\n", imax, acol.GetAt(0, imax))
		if imax > 0 {
			// pivot diagonal in symmetric matrix; will swap a11 and [imax,imax]
			applyPivotSym(&ABL, &ABR, imax, LOWER)
			p1.pivots[0] = imax + ATL.Rows() + 1
		} else {
			p1.pivots[0] = 0
		}
		if a11.Float() == 0.0 {
			err = onError("zero value on diagonal")
			return
		}

		//fmt.Printf("unblk pivoted %d, a11=%e, A:\n%v\n", imax, a11.Float(), A)
		//var Ablk matrix.FloatMatrix
		//merge1x2(&Ablk, &a21, &A22)
		//fmt.Printf("unblk update with a11=%e, a21|A22:\n%v\n", a11.Float(), &Ablk)

		// A22 = A22 - l21*d11*l21.T = A22 - a21*a21.T/a11; triangular update
		err = MVUpdateTrm(&A22, &a21, &a21, -1.0/a11.Float(), LOWER)

		// l21 = a21/a11
		InvScale(&a21, a11.Float())

		//merge1x2(&Ablk, &ABL, &ABR)
		//fmt.Printf("unblk imax=%d, Ablk:\n%v\n", imax, &Ablk)

		// ---------------------------------------------------------

		continue3x3to2x2(
			&ATL, &ATR,
			&ABL, &ABR, &A00, &a11, &A22, A, pBOTTOMRIGHT)
		contPivot3x1to2x1(
			&pT,
			&pB, &p0, &p1, p, pBOTTOM)
	}
	return
}

func findAndBuildPivot(AL, AR, WL, WR *matrix.FloatMatrix, k int) int {
	var dg, acol, wcol, wrow matrix.FloatMatrix

	// updated diagonal values on last column of workspace
	WR.SubMatrix(&dg, 0, WR.Cols()-1, AR.Rows(), 1)

	// find on-diagonal maximun value
	dmax := IAMax(&dg)
	//fmt.Printf("dmax=%d, val=%e\n", dmax, dg.GetAt(dmax, 0))

	// copy to first column of WR and update with factorized columns
	WR.SubMatrix(&wcol, 0, 0, WR.Rows(), 1)
	if dmax == 0 {
		AR.SubMatrix(&acol, 0, 0, AR.Rows(), 1)
		acol.CopyTo(&wcol)
	} else {
		AR.SubMatrix(&acol, dmax, 0, 1, dmax+1)
		acol.CopyTo(&wcol)
		if dmax < AR.Rows()-1 {
			var wrst matrix.FloatMatrix
			WR.SubMatrix(&wrst, dmax, 0, wcol.Rows()-dmax, 1)
			AR.SubMatrix(&acol, dmax, dmax, AR.Rows()-dmax, 1)
			acol.CopyTo(&wrst)
		}
	}
	if k > 0 {
		WL.SubMatrix(&wrow, dmax, 0, 1, WL.Cols())
		//fmt.Printf("update with wrow:%v\n", &wrow)
		//fmt.Printf("update wcol\n%v\n", &wcol)
		MVMult(&wcol, AL, &wrow, -1.0, 1.0, NOTRANS)
		//fmt.Printf("updated wcol:\n%v\n", &wcol)
	}
	if dmax > 0 {
		// pivot column in workspace
		t0 := WR.GetAt(0, 0)
		WR.SetAt(0, 0, WR.GetAt(dmax, 0))
		WR.SetAt(dmax, 0, t0)
		// pivot on diagonal
		t0 = dg.GetAt(0, 0)
		dg.SetAt(0, 0, dg.GetAt(dmax, 0))
		dg.SetAt(dmax, 0, t0)
	}
	return dmax
}

func unblkBoundedLowerLDL(A, W *matrix.FloatMatrix, p *pPivots, ncol int) (error, int) {
	var ATL, ATR, ABL, ABR matrix.FloatMatrix
	var A00, a10, a11, A20, a21, A22, adiag, wcol matrix.FloatMatrix
	var w00, w10, w11 matrix.FloatMatrix
	var pT, pB, p0, p1, p2 pPivots
	var err error = nil

	partition2x2(
		&ATL, &ATR,
		&ABL, &ABR, A, 0, 0, pTOPLEFT)
	partitionPivot2x1(
		&pT,
		&pB, p, 0, pTOP)

	// copy current diagonal to last column of workspace
	W.SubMatrix(&wcol, 0, W.Cols()-1, A.Rows(), 1)
	A.Diag(&adiag)
	adiag.CopyTo(&wcol)
	//fmt.Printf("initial diagonal:\n%v\n", &wcol)

	nc := 0
	for ABR.Cols() > 0 && nc < ncol {

		partition2x2(
			&w00, nil,
			&w10, &w11, W, nc, nc, pTOPLEFT)

		dmax := findAndBuildPivot(&ABL, &ABR, &w10, &w11, nc)
		//fmt.Printf("dmax=%d\n", dmax)
		if dmax > 0 {
			// pivot diagonal in symmetric matrix; will swap a11 [0,0] and [imax,imax]
			applyPivotSym(&ABL, &ABR, dmax, LOWER)
			swapRows(&w10, 0, dmax)
			pB.pivots[0] = dmax + ATL.Rows() + 1
		} else {
			pB.pivots[0] = 0
		}

		//fmt.Printf("blk pivoted %d, A:\n%v\nW:\n%v\n", dmax, A, W)
		repartition2x2to3x3(&ATL,
			&A00, nil, nil,
			&a10, &a11, nil,
			&A20, &a21, &A22, A, 1, pBOTTOMRIGHT)
		repartPivot2x1to3x1(&pT,
			&p0, &p1, &p2 /**/, p, 1, pBOTTOM)

		// --------------------------------------------------------

		// Copy updated column from working space
		w11.SubMatrix(&wcol, 1, 0, a21.Rows(), 1)
		wcol.CopyTo(&a21)
		a11.SetAt(0, 0, w11.GetAt(0, 0))
		// l21 = a21/a11
		InvScale(&a21, a11.Float())
		// here: wcol == l21*d11 == a21
		if ncol-nc > 1 {
			// update diagonal in workspace if not last column of block
			w11.SubMatrix(&adiag, 1, w11.Cols()-1, a21.Rows(), 1)
			MVUpdateDiag(&adiag, &wcol, &wcol, -1.0/a11.Float())
		}
		//fmt.Printf("nc=%d, a11=%e\n", nc, a11.Float())
		//fmt.Printf("l21\n%v\n", &a21)
		//fmt.Printf("a21\n%v\n", &wcol)
		//fmt.Printf("diag\n%v\n", &adiag)
		//var Ablk, wblk matrix.FloatMatrix
		//merge1x2(&Ablk, &ABL, &ABR)
		//merge1x2(&wblk, &w10, &w11)
		//fmt.Printf("unblk Ablk:\n%v\n", &Ablk)
		//fmt.Printf("unblk wblk:\n%v\n", &wblk)

		// ---------------------------------------------------------

		nc++
		continue3x3to2x2(
			&ATL, &ATR,
			&ABL, &ABR, &A00, &a11, &A22, A, pBOTTOMRIGHT)
		contPivot3x1to2x1(
			&pT,
			&pB, &p0, &p1, p, pBOTTOM)
	}
	return err, nc
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
func blkLowerLDL(A, W *matrix.FloatMatrix, p *pPivots, nb int) (err error) {
	var ATL, ATR, ABL, ABR matrix.FloatMatrix
	var A00, A10, A11, A20, A21, A22 matrix.FloatMatrix
	var /*D1,*/ wrk matrix.FloatMatrix
	var pT, pB, p0, p1, p2 pPivots
	var nblk int

	err = nil
	partition2x2(
		&ATL, &ATR,
		&ABL, &ABR, A, 0, 0, pTOPLEFT)
	partitionPivot2x1(
		&pT,
		&pB, p, 0, pTOP)

	for ABR.Cols() > nb {
		err, nblk = unblkBoundedLowerLDL(&ABR, W, &pB, nb)

		repartition2x2to3x3(&ATL,
			&A00, nil, nil,
			&A10, &A11, nil,
			&A20, &A21, &A22, A, nblk, pBOTTOMRIGHT)
		repartPivot2x1to3x1(&pT,
			&p0, &p1, &p2 /**/, p, nblk, pBOTTOM)

		// --------------------------------------------------------

		// wrk = D1*L21.T
		W.SubMatrix(&wrk, nblk, 0, A21.Rows(), nblk)

		// A22 = A22 - L21*D1*L21.T = A22 - L21*wrk.T
		UpdateTrm(&A22, &A21, &wrk, -1.0, 1.0, LOWER|TRANSB)

		applyRowPivots(&ABL, &p1, 0, FORWARD)
		scalePivots(&p1, ATL.Rows())

		//var Ablk matrix.FloatMatrix
		//merge1x2(&Ablk, &ABL, &ABR)
		//fmt.Printf("blk nblk=%d, Ablk:\n%v\n", nblk, &Ablk)

		// ---------------------------------------------------------

		continue3x3to2x2(
			&ATL, &ATR,
			&ABL, &ABR, &A00, &A11, &A22, A, pBOTTOMRIGHT)
		contPivot3x1to2x1(
			&pT,
			&pB, &p0, &p1, p, pBOTTOM)
	}
	if ABR.Cols() > 0 {
		// A11 = LDL(A11)
		err = unblkLowerLDL(&ABR, &pB)
		if err != nil {
			return
		}
		//fmt.Printf("unblk pivots: %v\n", pB.pivots)
		applyRowPivots(&ABL, &pB, 0, FORWARD)
		scalePivots(&pB, ATL.Rows())
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
func unblkUpperLDL(A *matrix.FloatMatrix, p *pPivots) (err error) {
	var ATL, ATR, ABL, ABR matrix.FloatMatrix
	var A00, a01, A02, a11, a12, A22 matrix.FloatMatrix
	var AL, AR, acol matrix.FloatMatrix
	var pT, pB, p0, p1, p2 pPivots

	err = nil
	partition2x2(
		&ATL, &ATR,
		&ABL, &ABR, A, 0, 0, pBOTTOMRIGHT)
	partitionPivot2x1(
		&pT,
		&pB, p, 0, pBOTTOM)

	for ATL.Rows() > 0 {
		repartition2x2to3x3(&ATL,
			&A00, &a01, &A02,
			nil, &a11, &a12,
			nil, nil, &A22, A, 1, pTOPLEFT)
		repartPivot2x1to3x1(&pT,
			&p0, &p1, &p2 /**/, p, 1, pTOP)

		// --------------------------------------------------------
		// search diagonal; diag(A00;a11)
		ATL.Diag(&acol)
		//merge2x1(&acol, &a01, &a11)
		imax := IAMax(&acol)
		if imax < ATL.Rows()-1 {
			merge1x2(&AL, &ATL, &ATR)
			merge1x2(&AR, &a11, &a12)
			// pivot diagonal in symmetric matrix; will swap a11 and [imax,imax]
			applyPivotSym(&AL, &AR, imax, UPPER)
			p1.pivots[0] = imax + 1
		} else {
			p1.pivots[0] = 0
		}

		if a11.Float() == 0.0 {
			err = onError("zero on diagonal.")
			return
		}
		// A00 = A00 - u01*d11*u01.T = A00 - a01*a01.T/a11; triangular update
		err = MVUpdateTrm(&A00, &a01, &a01, -1.0/a11.Float(), UPPER)

		// u01 = a01/a11
		InvScale(&a01, a11.Float())
		// ---------------------------------------------------------

		continue3x3to2x2(
			&ATL, &ATR,
			&ABL, &ABR, &A00, &a11, &A22, A, pTOPLEFT)
		contPivot3x1to2x1(
			&pT,
			&pB, &p0, &p1, p, pTOP)
	}
	return
}

func blkUpperLDL(A, W *matrix.FloatMatrix, p *pPivots, nb int) (err error) {
	var ATL, ATR, ABL, ABR matrix.FloatMatrix
	var A00, A01, A02, A11, A12, A22 matrix.FloatMatrix
	var D1, wrk matrix.FloatMatrix
	var pT, pB, p0, p1, p2 pPivots

	err = nil
	partition2x2(
		&ATL, &ATR,
		&ABL, &ABR, A, 0, 0, pBOTTOMRIGHT)
	partitionPivot2x1(
		&pT,
		&pB, p, 0, pBOTTOM)

	for ATL.Rows() > 0 {
		repartition2x2to3x3(&ATL,
			&A00, &A01, &A02,
			nil, &A11, &A12,
			nil, nil, &A22, A, nb, pTOPLEFT)
		repartPivot2x1to3x1(&pT,
			&p0, &p1, &p2 /**/, p, nb, pTOP)

		// --------------------------------------------------------

		// A11 = LDL(A11)
		err = unblkUpperLDL(&A11, &p1)
		if err != nil {
			return
		}
		applyColPivots(&A01, &p1, 0, BACKWARD)
		applyRowPivots(&A12, &p1, 0, BACKWARD)
		scalePivots(&p1, ATL.Rows()-A11.Rows())

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
		contPivot3x1to2x1(
			&pT,
			&pB, &p0, &p1, p, pTOP)
	}
	return
}

/*
 * Compute an LDLT factorization of a symmetric N-by-N matrix with partial pivoting.
 *
 * Arguments:
 *   A      On entry, the N-by-N matrix to be factored. On exit the factor
 *          L and 1-by-1 diagonal D from factorization A = L*D*L.T, the unit diagonal
 *          of L are not stored. Or the factor U and diagonal D from factorization
 *          A = U*D*U.T if flag bit UPPER is set.
 *
 *   W      Work space for blocking invocations, matrix of size N-by-nb.
 *
 *   ipiv   Pivot indeces, for each non-zero element ipiv[k] the k'th row is exchanged with
 *          ipiv[k]-1'th row.
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
func DecomposeLDL(A, W *matrix.FloatMatrix, ipiv []int, flags Flags, nb int) (*matrix.FloatMatrix, error) {
	var err error
	if A.Cols() != A.Rows() {
		return nil, errors.New("A not a square matrix")
	}
	for k, _ := range ipiv {
		ipiv[k] = 0
	}
	if A.Cols() < nb || nb == 0 {
		if flags&LOWER != 0 {
			err = unblkLowerLDL(A, &pPivots{ipiv})
		} else {
			err = unblkUpperLDL(A, &pPivots{ipiv})
		}
	} else {
		if flags&LOWER != 0 {
			err = blkLowerLDL(A, W, &pPivots{ipiv}, nb)
		} else {
			err = blkUpperLDL(A, W, &pPivots{ipiv}, nb)
		}
	}
	return A, err
}

/*
 * Solves a system system of linear equations A*X = B with symmetric positive
 * definite matrix A using the LDL factorization A = L*D*L.T or A = U*D*U.T
 * computed by DecomposeLDL().
 *
 * Arguments:
 *  B      On entry, the unpermuted right hand side matrix B.
 *         On exit, the solution matrix X.
 *
 *  A      The triangular factor U or L from LDL factorization as computed by
 *         DecomposeLDL().
 *
 *  ipiv   Pivot indeces, for each non-zero element ipiv[k] the k'th row is
 *         exchanged with ipiv[k]-1'th row.
 *
 *  flags  Indicator of which factor is stored in A. If flags&UPPER then upper
 *         triangle of A is stored. If flags&LOWER then lower triangle of A is
 *         stored.
 *
 * Notes:
 *  On entry matrix B is permuted according ipiv vector and on exit
 *  rearraged to original row order.
 */
func SolveLDL(B, A *matrix.FloatMatrix, ipiv []int, flags Flags) {
	if flags&UPPER != 0 {
		// X = (U*D*U.T).-1*B => U.-T*(D.-1*(U.-1*B))
		// arrange to match factorization
		applyRowPivots(B, &pPivots{ipiv}, 0, BACKWARD)
		// solve
		SolveTrm(B, A, 1.0, UPPER|UNIT)
		SolveDiag(B, A, LEFT)
		SolveTrm(B, A, 1.0, UPPER|UNIT|TRANSA)
		// rearrange to original
		applyRowPivots(B, &pPivots{ipiv}, 0, FORWARD)

	} else if flags&LOWER != 0 {
		// X = (L*D*L.T).-1*B = L.-T*(D*-1(L.-1*B))
		// arrange to match factorization
		applyRowPivots(B, &pPivots{ipiv}, 0, FORWARD)
		// solve
		SolveTrm(B, A, 1.0, LOWER|UNIT)
		SolveDiag(B, A, LEFT)
		SolveTrm(B, A, 1.0, LOWER|UNIT|TRANSA)
		// rearrange to original
		applyRowPivots(B, &pPivots{ipiv}, 0, BACKWARD)
	}
}

// Local Variables:
// tab-width: 4
// indent-tabs-mode: nil
// End:
