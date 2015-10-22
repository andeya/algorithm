// Copyright (c) Harri Rautila, 2013

// This file is part of github.com/hrautila/matops package. It is free software,
// distributed under the terms of GNU Lesser General Public License Version 3, or
// any later version. See the COPYING tile included in this archive.

package matops

import (
	"errors"
	"github.com/henrylee2cn/algorithm/matrix"
	"math"
	//"fmt"
)

const bkALPHA = 0.6403882032022075 // (1.0 + sqrt(17.0))/8.0

/*
 * Apply diagonal pivot (row and column swapped) to symmetric matrix blocks.
 *
 * LOWER triangular; moving from top-left to bottom-right
 *
 *    -----------------------
 *    | d
 *    | x P1 x  x  x  P2     -- current row/col 'srcix'
 *    | x S2 d  x  x  x
 *    | x S2 x  d  x  x
 *    | x S2 x  x  d  x
 *    | x P2 D2 D2 D2 P3     -- swap with row/col 'dstix'
 *    | x S3 x  x  x  D3 d
 *    | x S3 x  x  x  D3 x d
 *         (AR)
 *
 * UPPER triangular; moving from bottom-right to top-left
 *
 *    d x D3 x  x  x  S3 x |
 *      d D3 x  x  x  S3 x |
 *        P3 D2 D2 D2 P2 x |  -- dstinx
 *           d  x  x  S2 x |
 *              d  x  S2 x |
 *                 d  S2 x |
 *                    P1 x |  -- srcinx
 *                       d |
 *    ----------------------
 *               (ABR)
 */
func applyBKPivotSym(AR *matrix.FloatMatrix, srcix, dstix int, flags Flags) {
	var s, d matrix.FloatMatrix
	if flags&LOWER != 0 {
		// S2 -- D2
		AR.SubMatrix(&s, srcix+1, srcix, dstix-srcix-1, 1)
		AR.SubMatrix(&d, dstix, srcix+1, 1, dstix-srcix-1)
		Swap(&s, &d)
		// S3 -- D3
		AR.SubMatrix(&s, dstix+1, srcix, AR.Rows()-dstix-1, 1)
		AR.SubMatrix(&d, dstix+1, dstix, AR.Rows()-dstix-1, 1)
		Swap(&s, &d)
		// swap P1 and P3
		p1 := AR.GetAt(srcix, srcix)
		p3 := AR.GetAt(dstix, dstix)
		AR.SetAt(srcix, srcix, p3)
		AR.SetAt(dstix, dstix, p1)
		return
	}
	if flags&UPPER != 0 {
		// AL is ATL, AR is ATR; P1 is AL[srcix, srcix];
		// S2 -- D2
		AR.SubMatrix(&s, dstix+1, srcix, srcix-dstix-1, 1)
		AR.SubMatrix(&d, dstix, dstix+1, 1, srcix-dstix-1)
		Swap(&s, &d)
		// S3 -- D3
		AR.SubMatrix(&s, 0, srcix, dstix, 1)
		AR.SubMatrix(&d, 0, dstix, dstix, 1)
		Swap(&s, &d)
		//fmt.Printf("3, AR=%v\n", AR)
		// swap P1 and P3
		p1 := AR.GetAt(srcix, srcix)
		p3 := AR.GetAt(dstix, dstix)
		AR.SetAt(srcix, srcix, p3)
		AR.SetAt(dstix, dstix, p1)
		return
	}
}

/*
 * α = (1 + sqrt(17))/8
 * λ = |a(r,1)| = max{|a(2,1)|, . . . , |a(m,1)|}
 * if λ > 0
 *     if |a(1,1)| ≥ αλ
 *         use a11 as 1-by-1 pivot
 *     else
 *         σ = |a(p,r)| = max{|a(1,r)|,..., |a(r−1,r)|, |a(r+1,r)|,..., |a(m,r)|}
 *         if |a(1,1) |σ ≥ αλ^2
 *             use a(1,1) as 1-by-1 pivot
 *         else if |a(r,r)| ≥ ασ
 *             use a(r,r) as 1-by-1 pivot
 *         else
 *                  a11 | ar1
 *             use  --------  as 2-by-2 pivot
 *                  ar1 | arr
 *         end
 *     end
 * end
 */

func findBKPivot(A *matrix.FloatMatrix, flags Flags) (int, int) {
	var r, q int
	var rcol, qrow matrix.FloatMatrix
	if flags&LOWER != 0 {
		if A.Rows() == 1 {
			return 0, 1
		}
		amax := math.Abs(A.GetAt(0, 0))
		// column below diagonal at [0, 0]
		A.SubMatrix(&rcol, 1, 0, A.Rows()-1, 1)
		r = IAMax(&rcol) + 1
		// max off-diagonal on first column at index r
		rmax := math.Abs(A.GetAt(r, 0))
		//fmt.Printf("m(A)=%d, r=%d, rmax=%e, amax=%e\n", m(A), r, rmax, amax)
		if amax >= bkALPHA*rmax {
			// no pivoting, 1x1 diagonal
			return 0, 1
		}
		// max off-diagonal on r'th row at index q
		A.SubMatrix(&qrow, r, 0, 1, r /*+1*/)
		q = IAMax(&qrow)
		qmax := math.Abs(A.GetAt(r, q /*+1*/))
		if r < A.Rows()-1 {
			// rest of the r'th row after diagonal
			A.SubMatrix(&qrow, r+1, r, A.Rows()-r-1, 1)
			q = IAMax(&qrow)
			//fmt.Printf("qrow: %d, q: %d\n", qrow.NumElements(), q)
			qmax2 := math.Abs(qrow.GetAt(q, 0))
			if qmax2 > qmax {
				qmax = qmax2
			}
		}
		//fmt.Printf("m(A)=%d: q=%d, qmax=%e %v\n", m(A), q, qmax, &qrow)
		//arr := math.Abs(A.GetAt(r, r))
		//fmt.Printf("unblk: r=%d, q=%d, amax=%e, rmax=%e, qmax=%e, Arr=%e\n", r, q, amax, rmax, qmax, arr)

		if amax >= bkALPHA*rmax*(rmax/qmax) {
			// no pivoting, 1x1 diagonal
			return 0, 1
		}
		if math.Abs(A.GetAt(r, r)) >= bkALPHA*qmax {
			// 1x1 pivoting and interchange with k, r
			return r, 1
		} else {
			// 2x2 pivoting and interchange with k+1, r
			return r, 2
		}
	}
	if flags&UPPER != 0 {
		if A.Rows() == 1 {
			return 0, 1
		}
		//fmt.Printf("upper A:\n%v\n", A)
		lastcol := A.Rows() - 1
		amax := math.Abs(A.GetAt(lastcol, lastcol))
		// column above [A.Rows()-1, A.Rows()-1]
		A.SubMatrix(&rcol, 0, lastcol, lastcol, 1)
		r = IAMax(&rcol)
		// max off-diagonal on first column at index r
		rmax := math.Abs(A.GetAt(r, lastcol))
		//fmt.Printf("m(A)=%d, r=%d, rmax=%e, amax=%e\n", m(A), r, rmax, amax)
		if amax >= bkALPHA*rmax {
			// no pivoting, 1x1 diagonal
			return -1, 1
		}
		// max off-diagonal on r'th row at index q
		//  a) rest of the r'th row above diagonal
		qmax := 0.0
		if r > 0 {
			A.SubMatrix(&qrow, 0, r, r, 1)
			q = IAMax(&qrow)
			qmax = math.Abs(A.GetAt(q, r /*+1*/))
		}
		//  b) elements right of diagonal
		A.SubMatrix(&qrow, r, r+1, 1, lastcol-r)
		q = IAMax(&qrow)
		//fmt.Printf("qrow: %d, q: %d, data: %v\n", qrow.NumElements(), q, &qrow)
		qmax2 := math.Abs(qrow.GetAt(0, q))
		if qmax2 > qmax {
			qmax = qmax2
		}

		//fmt.Printf("m(A)=%d: q=%d, qmax=%e %v\n", m(A), q, qmax, &qrow)
		//fmt.Printf("unblk: r=%d, q=%d, amax=%e, rmax=%e, qmax=%e\n", r, q, amax, rmax, qmax)

		if amax >= bkALPHA*rmax*(rmax/qmax) {
			// no pivoting, 1x1 diagonal
			return -1, 1
		}
		if math.Abs(A.GetAt(r, r)) >= bkALPHA*qmax {
			// 1x1 pivoting and interchange with k, r
			return r, 1
		} else {
			// 2x2 pivoting and interchange with k+1, r
			return r, 2
		}
	}
	return 0, 1
}

/*
 * Unblocked Bunch-Kauffman LDL factorization.
 *
 * Corresponds lapack.DSYTF2
 */
func unblkDecompBKLower(A, wrk *matrix.FloatMatrix, p *pPivots) (error, int) {
	var err error
	var ATL, ATR, ABL, ABR matrix.FloatMatrix
	var A00, a10t, a11, A20, a21, A22, a11inv matrix.FloatMatrix
	var cwrk matrix.FloatMatrix
	var pT, pB, p0, p1, p2 pPivots

	err = nil
	nc := 0

	partition2x2(
		&ATL, &ATR,
		&ABL, &ABR, A, 0, 0, pTOPLEFT)
	partitionPivot2x1(
		&pT,
		&pB, p, 0, pTOP)

	// permanent working space for symmetric inverse of a11
	wrk.SubMatrix(&a11inv, 0, wrk.Cols()-2, 2, 2)
	a11inv.SetAt(1, 0, -1.0)
	a11inv.SetAt(0, 1, -1.0)

	for ABR.Cols() > 0 {

		r, np := findBKPivot(&ABR, LOWER)
		if r != 0 && r != np-1 {
			// pivoting needed; do swaping here
			applyBKPivotSym(&ABR, np-1, r, LOWER)
			if np == 2 {
				/*
				 *          [0,0] | [r,0]
				 * a11 ==   -------------  2-by-2 pivot, swapping [1,0] and [r,0]
				 *          [r,0] | [r,r]
				 */
				t := ABR.GetAt(1, 0)
				ABR.SetAt(1, 0, ABR.GetAt(r, 0))
				ABR.SetAt(r, 0, t)
			}
			//fmt.Printf("unblk: ABR after %d pivot [r=%d]:\n%v\n", np, r, &ABR)
		}

		// repartition according the pivot size
		repartition2x2to3x3(&ATL,
			&A00, nil, nil,
			&a10t, &a11, nil,
			&A20, &a21, &A22 /**/, A, np, pBOTTOMRIGHT)
		repartPivot2x1to3x1(&pT,
			&p0,
			&p1,
			&p2 /**/, p, np, pBOTTOM)
		// ------------------------------------------------------------

		if np == 1 {
			// A22 = A22 - a21*a21.T/a11
			MVUpdateTrm(&A22, &a21, &a21, -1.0/a11.Float(), LOWER)
			// a21 = a21/a11
			InvScale(&a21, a11.Float())
			// store pivot point relative to original matrix
			p1.pivots[0] = r + ATL.Rows() + 1
		} else if np == 2 {
			/* from Bunch-Kaufmann 1977:
			 *  (E2 C.T) = ( I2      0      )( E  0      )( I[n-2] E.-1*C.T )
			 *  (C  B  )   ( C*E.-1  I[n-2] )( 0  A[n-2] )( 0      I2       )
			 *
			 *  A[n-2] = B - C*E.-1*C.T
			 *
			 *  E.-1 is inverse of a symmetric matrix, cannot use
			 *  triangular solve. We calculate inverse of 2x2 matrix.
			 *  Following is inspired by lapack.SYTF2
			 *
			 *      a | b      1        d | -b         b         d/b | -1
			 *  inv ----- =  ------  * ------  =  ----------- * --------
			 *      b | d    (ad-b^2)  -b |  a    (a*d - b^2)     -1 | a/b
			 *
			 */
			a := a11.GetAt(0, 0)
			b := a11.GetAt(1, 0)
			d := a11.GetAt(1, 1)
			a11inv.SetAt(0, 0, d/b)
			a11inv.SetAt(1, 1, a/b)
			// denominator: (a/b)*(d/b)-1.0 == (a*d - b^2)/b^2
			scale := 1.0 / ((a/b)*(d/b) - 1.0)
			scale /= b

			// cwrk = a21
			wrk.SubMatrix(&cwrk, 2, 0, a21.Rows(), a21.Cols())
			a21.CopyTo(&cwrk)
			// a21 = a21*a11.-1
			Mult(&a21, &cwrk, &a11inv, scale, 0.0, NOTRANS)
			// A22 = A22 - a21*a11.-1*a21.T = A22 - a21*cwrk.T
			UpdateTrm(&A22, &a21, &cwrk, -1.0, 1.0, LOWER|TRANSB)

			// store pivot point relative to original matrix
			p1.pivots[0] = -(r + ATL.Rows() + 1)
			p1.pivots[1] = p1.pivots[0]
		}

		/*
		   if m(&ABR) < 5 {
		       var Ablk matrix.FloatMatrix
		       merge1x2(&Ablk, &ABL, &ABR)
		       fmt.Printf("unblocked EOL: Ablk r=%d, nc=%d. np=%d\n%v\n", r, nc, np, &Ablk)
		   }
		*/
		// ------------------------------------------------------------
		nc += np
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
 * Unblocked Bunch-Kauffman LDL factorization.
 *
 * Corresponds lapack.DSYTF2
 */
func unblkDecompBKUpper(A, wrk *matrix.FloatMatrix, p *pPivots) (error, int) {
	var err error
	var ATL, ATR, ABL, ABR matrix.FloatMatrix
	var A00, a01, A02, a12t, a11, A22, a11inv matrix.FloatMatrix
	var cwrk matrix.FloatMatrix
	var pT, pB, p0, p1, p2 pPivots

	err = nil
	nc := 0

	partition2x2(
		&ATL, &ATR,
		&ABL, &ABR, A, 0, 0, pBOTTOMRIGHT)
	partitionPivot2x1(
		&pT,
		&pB, p, 0, pBOTTOM)

	// permanent working space for symmetric inverse of a11
	wrk.SubMatrix(&a11inv, 0, wrk.Cols()-2, 2, 2)
	a11inv.SetAt(1, 0, -1.0)
	a11inv.SetAt(0, 1, -1.0)

	for ATL.Cols() > 0 {

		nr := ATL.Rows() - 1
		r, np := findBKPivot(&ATL, UPPER)
		if r != -1 /*&& r != np-1*/ {
			// pivoting needed; do swaping here
			//fmt.Printf("pre-pivot ATL [%d]:\n%v\n", ATL.Rows()-np, &ATL)
			applyBKPivotSym(&ATL, ATL.Rows()-np, r, UPPER)
			if np == 2 {
				/*
				 *         [r,r] | [r, nr]
				 * a11 ==  ---------------  2-by-2 pivot, swapping [nr-1,nr] and [r,nr]
				 *         [r,0] | [nr,nr]
				 */
				t := ATL.GetAt(nr-1, nr)
				ATL.SetAt(nr-1, nr, ATL.GetAt(r, nr))
				ATL.SetAt(r, nr, t)
			}
			//fmt.Printf("unblk: ATL after %d pivot [r=%d]:\n%v\n", np, r, &ATL)
		}

		// repartition according the pivot size
		repartition2x2to3x3(&ATL,
			&A00, &a01, &A02,
			nil, &a11, &a12t,
			nil, nil, &A22 /**/, A, np, pTOPLEFT)
		repartPivot2x1to3x1(&pT,
			&p0,
			&p1,
			&p2 /**/, p, np, pTOP)
		// ------------------------------------------------------------

		if np == 1 {
			// A00 = A00 - a01*a01.T/a11
			MVUpdateTrm(&A00, &a01, &a01, -1.0/a11.Float(), UPPER)
			// a01 = a01/a11
			InvScale(&a01, a11.Float())
			if r == -1 {
				p1.pivots[0] = ATL.Rows()
			} else {
				p1.pivots[0] = r + 1
			}
		} else if np == 2 {
			/*
			 * See comments on unblkDecompBKLower().
			 */
			a := a11.GetAt(0, 0)
			b := a11.GetAt(0, 1)
			d := a11.GetAt(1, 1)
			a11inv.SetAt(0, 0, d/b)
			a11inv.SetAt(1, 1, a/b)
			// denominator: (a/b)*(d/b)-1.0 == (a*d - b^2)/b^2
			scale := 1.0 / ((a/b)*(d/b) - 1.0)
			scale /= b

			// cwrk = a21
			wrk.SubMatrix(&cwrk, 2, 0, a01.Rows(), a01.Cols())
			a01.CopyTo(&cwrk)
			//fmt.Printf("cwrk:\n%v\n", &cwrk)
			//fmt.Printf("a11inv:\n%v\n", &a11inv)
			// a01 = a01*a11.-1
			Mult(&a01, &cwrk, &a11inv, scale, 0.0, NOTRANS)
			// A00 = A00 - a01*a11.-1*a01.T = A00 - a01*cwrk.T
			UpdateTrm(&A00, &a01, &cwrk, -1.0, 1.0, UPPER|TRANSB)

			p1.pivots[0] = -(r + 1)
			p1.pivots[1] = p1.pivots[0]
		}

		// ------------------------------------------------------------
		nc += np
		continue3x3to2x2(
			&ATL, &ATR,
			&ABL, &ABR, &A00, &a11, &A22, A, pTOPLEFT)
		contPivot3x1to2x1(
			&pT,
			&pB, &p0, &p1, p, pTOP)

	}
	return err, nc
}

/*
 * Find diagonal pivot and build incrementaly updated block.
 *
 *  (AL)  (AR)                   (WL)  (WR)
 *  --------------------------   ----------    k'th row in W
 *  x x | c1                     w w | k kp1
 *  x x | c1 d                   w w | k kp1
 *  x x | c1 x  d                w w | k kp1
 *  x x | c1 x  x  d             w w | k kp1
 *  x x | c1 r2 r2 r2 r2         w w | k kp1
 *  x x | c1 x  x  x  r2 d       w w | k kp1
 *  x x | c1 x  x  x  r2 x d     w w | k kp1
 *
 * Matrix AR contains the unfactored part of the matrix and AL the already
 * factored columns. Matrix WL is updated values of factored part ie.
 * w(i) = l(i)d(i). Matrix WR will have updated values for next column.
 * Column WR(k) contains updated AR(c1) and WR(kp1) possible pivot row AR(r2).
 *
 *
 */
func findAndBuildBKPivotLower(AL, AR, WL, WR *matrix.FloatMatrix, k int) (int, int) {
	var r, q int
	var rcol, qrow, src, wk, wkp1, wrow matrix.FloatMatrix

	// Copy AR column 0 to WR column 0 and update with WL[0:]
	AR.SubMatrix(&src, 0, 0, AR.Rows(), 1)
	WR.SubMatrix(&wk, 0, 0, AR.Rows(), 1)
	src.CopyTo(&wk)
	if k > 0 {
		WL.SubMatrix(&wrow, 0, 0, 1, WL.Cols())
		MVMult(&wk, AL, &wrow, -1.0, 1.0, NOTRANS)
		//fmt.Printf("wk after update:\n%v\n", &wk)
	}
	if AR.Rows() == 1 {
		return 0, 1
	}
	amax := math.Abs(WR.GetAt(0, 0))

	// find max off-diagonal on first column.
	WR.SubMatrix(&rcol, 1, 0, AR.Rows()-1, 1)
	//fmt.Printf("rcol:\n%v\n", &rcol)
	// r is row index and rmax is its absolute value
	r = IAMax(&rcol) + 1
	rmax := math.Abs(rcol.GetAt(r-1, 0))
	//fmt.Printf("r=%d, amax=%e, rmax=%e\n", r, amax, rmax)
	if amax >= bkALPHA*rmax {
		// no pivoting, 1x1 diagonal
		return 0, 1
	}
	// Now we need to copy row r to WR[:,1] and update it
	WR.SubMatrix(&wkp1, 0, 1, AR.Rows(), 1)
	AR.SubMatrix(&qrow, r, 0, 1, r+1)
	qrow.CopyTo(&wkp1)
	//fmt.Printf("m(AR)=%d, r=%d, qrow: %v\n", AR.Rows(), r, &qrow)
	if r < AR.Rows()-1 {
		var wkr matrix.FloatMatrix
		AR.SubMatrix(&qrow, r, r, AR.Rows()-r, 1)
		wkp1.SubMatrix(&wkr, r, 0, wkp1.Rows()-r, 1)
		qrow.CopyTo(&wkr)
		//fmt.Printf("m(AR)=%d, r=%d, qrow: %v\n", AR.Rows(), r, &qrow)
	}
	if k > 0 {
		// update wkp1
		WL.SubMatrix(&wrow, r, 0, 1, WL.Cols())
		//fmt.Printf("initial wpk1:\n%v\n", &wkp1)
		MVMult(&wkp1, AL, &wrow, -1.0, 1.0, NOTRANS)
		//fmt.Printf("updated wpk1:\n%v\n", &wkp1)
	}

	// set on-diagonal entry to zero to avoid finding it
	p1 := wkp1.GetAt(r, 0)
	wkp1.SetAt(r, 0, 0.0)
	// max off-diagonal on r'th column/row at index q
	q = IAMax(&wkp1)
	qmax := math.Abs(wkp1.GetAt(q, 0))
	// restore on-diagonal entry
	wkp1.SetAt(r, 0, p1)
	//arr := math.Abs(WR.GetAt(r, 1))
	//fmt.Printf("blk: r=%d, q=%d, amax=%e, rmax=%e, qmax=%e, Arr=%e\n", r, q, amax, rmax, qmax, arr)

	if amax >= bkALPHA*rmax*(rmax/qmax) {
		// no pivoting, 1x1 diagonal
		return 0, 1
	}
	// if q == r then qmax is not off-diagonal, qmax == WR[r,1] and
	// we get 1x1 pivot as following is always true
	if math.Abs(WR.GetAt(r, 1)) >= bkALPHA*qmax {
		// 1x1 pivoting and interchange with k, r
		// pivot row in column WR[:,1] to W[:,0]
		//pr := WR.GetAt(r, 1)
		//_ = pr
		WR.SubMatrix(&src, 0, 1, AR.Rows(), 1)
		WR.SubMatrix(&wkp1, 0, 0, AR.Rows(), 1)
		src.CopyTo(&wkp1)
		wkp1.SetAt(0, 0, src.GetAt(r, 0))
		wkp1.SetAt(r, 0, src.GetAt(0, 0))
		return r, 1
	} else {
		// 2x2 pivoting and interchange with k+1, r
		return r, 2
	}
	return 0, 1
}

/*
 * Unblocked, bounded Bunch-Kauffman LDL factorization for at most ncol columns.
 * At most ncol columns are factorized and trailing matrix updates are restricted
 * to ncol columns. Also original columns are accumulated to working matrix, which
 * is used by calling blocked algorithm to update the trailing matrix with BLAS3
 * update.
 *
 * Corresponds lapack.DLASYF
 */
func unblkBoundedBKLower(A, wrk *matrix.FloatMatrix, p *pPivots, ncol int) (error, int) {
	var err error
	var ATL, ATR, ABL, ABR matrix.FloatMatrix
	var A00, a10t, a11, A20, a21, A22, a11inv matrix.FloatMatrix
	var w00, w10, w11 matrix.FloatMatrix
	var cwrk matrix.FloatMatrix
	//var s, d matrix.FloatMatrix
	var pT, pB, p0, p1, p2 pPivots

	err = nil
	nc := 0
	if ncol > A.Cols() {
		ncol = A.Cols()
	}

	partition2x2(
		&ATL, &ATR,
		&ABL, &ABR, A, 0, 0, pTOPLEFT)
	partitionPivot2x1(
		&pT,
		&pB, p, 0, pTOP)

	// permanent working space for symmetric inverse of a11
	wrk.SubMatrix(&a11inv, 0, wrk.Cols()-2, 2, 2)
	a11inv.SetAt(1, 0, -1.0)
	a11inv.SetAt(0, 1, -1.0)

	for ABR.Cols() > 0 && nc < ncol {

		partition2x2(
			&w00, nil,
			&w10, &w11, wrk, nc, nc, pTOPLEFT)

		//fmt.Printf("ABR:\n%v\n", &ABR)
		r, np := findAndBuildBKPivotLower(&ABL, &ABR, &w10, &w11, nc)
		//fmt.Printf("after find: r=%d, np=%d, ncol=%d, nc=%d\n", r, np, ncol, nc)
		if np > ncol-nc {
			// next pivot does not fit into ncol columns, restore last column,
			// return with number of factorized columns
			//fmt.Printf("np > ncol-nc: %d > %d\n", np, ncol-nc)
			return err, nc
			//goto undo
		}
		if r != 0 && r != np-1 {
			// pivoting needed; do swaping here
			applyBKPivotSym(&ABR, np-1, r, LOWER)
			// swap left hand rows to get correct updates
			swapRows(&ABL, np-1, r)
			swapRows(&w10, np-1, r)
			//ABL.SubMatrix(&s, np-1, 0, 1, ABL.Cols())
			//ABL.SubMatrix(&d, r,    0, 1, ABL.Cols())
			//Swap(&s, &d)
			//w10.SubMatrix(&s, np-1, 0, 1, w10.Cols())
			//w10.SubMatrix(&d, r,    0, 1, w10.Cols())
			//Swap(&s, &d)
			if np == 2 {
				/*
				 *          [0,0] | [r,0]
				 * a11 ==   -------------  2-by-2 pivot, swapping [1,0] and [r,0]
				 *          [r,0] | [r,r]
				 */
				t0 := w11.GetAt(1, 0)
				tr := w11.GetAt(r, 0)
				//fmt.Printf("nc=%d, t0=%e, tr=%e\n", nc, t0, tr)
				w11.SetAt(1, 0, tr)
				w11.SetAt(r, 0, t0)
				// interchange diagonal entries on w11[:,1]
				t0 = w11.GetAt(1, 1)
				tr = w11.GetAt(r, 1)
				w11.SetAt(1, 1, tr)
				w11.SetAt(r, 1, t0)
			}
			//fmt.Printf("pivoted A:\n%v\n", A)
			//fmt.Printf("pivoted wrk:\n%v\n", wrk)
		}

		// repartition according the pivot size
		repartition2x2to3x3(&ATL,
			&A00, nil, nil,
			&a10t, &a11, nil,
			&A20, &a21, &A22 /**/, A, np, pBOTTOMRIGHT)
		repartPivot2x1to3x1(&pT,
			&p0,
			&p1,
			&p2 /**/, p, np, pBOTTOM)
		// ------------------------------------------------------------

		if np == 1 {
			//
			w11.SubMatrix(&cwrk, np, 0, a21.Rows(), np)
			a11.SetAt(0, 0, w11.GetAt(0, 0))
			// a21 = a21/a11
			//fmt.Printf("np == 1: pre-update a21\n%v\n", &a21)
			cwrk.CopyTo(&a21)
			InvScale(&a21, a11.Float())
			//fmt.Printf("np == 1: cwrk\n%v\na21\n%v\n", &cwrk, &a21)
			// store pivot point relative to original matrix
			p1.pivots[0] = r + ATL.Rows() + 1
		} else if np == 2 {
			/*
			 * See comments for this block in unblkDecompBKLower().
			 */
			a := w11.GetAt(0, 0)
			b := w11.GetAt(1, 0)
			d := w11.GetAt(1, 1)
			a11inv.SetAt(0, 0, d/b)
			a11inv.SetAt(1, 1, a/b)
			// denominator: (a/b)*(d/b)-1.0 == (a*d - b^2)/b^2
			scale := 1.0 / ((a/b)*(d/b) - 1.0)
			scale /= b

			w11.SubMatrix(&cwrk, np, 0, a21.Rows(), np)
			// a21 = a21*a11.-1
			Mult(&a21, &cwrk, &a11inv, scale, 0.0, NOTRANS)
			a11.SetAt(0, 0, a)
			a11.SetAt(1, 0, b)
			a11.SetAt(1, 1, d)

			// store pivot point relative to original matrix
			p1.pivots[0] = -(r + ATL.Rows() + 1)
			p1.pivots[1] = p1.pivots[0]
		}

		/*
		   if m(&ABR) < 5 {
		       var Ablk, wblk, w5 matrix.FloatMatrix
		       merge1x2(&Ablk, &ABL, &ABR)
		       merge1x2(&wblk, &w10, &w11)
		       wblk.SubMatrix(&w5, 0, 0, Ablk.Rows(), wblk.Cols())
		       fmt.Printf("blocked EOL: Ablk r=%d, nc=%d. np=%d\n%v\n", r, nc, np, &Ablk)
		       fmt.Printf("wblk m(wblk)=%d:\n%v\n", m(&w5), &w5)
		   }
		*/
		// ------------------------------------------------------------
		nc += np
		continue3x3to2x2(
			&ATL, &ATR,
			&ABL, &ABR, &A00, &a11, &A22, A, pBOTTOMRIGHT)
		contPivot3x1to2x1(
			&pT,
			&pB, &p0, &p1, p, pBOTTOM)

	}
	// undo applied partial row pivots (AL, w00)
	//undo:
	return err, nc
}

func blkDecompBKLower(A, W *matrix.FloatMatrix, p *pPivots, nb int) (err error) {
	var ATL, ATR, ABL, ABR matrix.FloatMatrix
	var A00, A10, A11, A20, A21, A22 matrix.FloatMatrix
	var wrk matrix.FloatMatrix
	var pT, pB, p0, p1, p2 pPivots
	var nblk int = 0

	err = nil
	partition2x2(
		&ATL, &ATR,
		&ABL, &ABR, A, 0, 0, pTOPLEFT)
	partitionPivot2x1(
		&pT,
		&pB, p, 0, pTOP)

	for ABR.Cols() >= nb {
		err, nblk = unblkBoundedBKLower(&ABR, W, &pB, nb)

		// repartition nblk size
		repartition2x2to3x3(&ATL,
			&A00, nil, nil,
			&A10, &A11, nil,
			&A20, &A21, &A22, A, nblk, pBOTTOMRIGHT)
		repartPivot2x1to3x1(&pT,
			&p0, &p1, &p2 /**/, p, nblk, pBOTTOM)

		// --------------------------------------------------------
		// here [A11;A21] has been decomposed by unblkBoundedBKLower()
		// Now we need update A22

		// wrk is original A21
		W.SubMatrix(&wrk, nblk, 0, A21.Rows(), nblk)

		// A22 = A22 - L21*D1*L21.T = A22 - L21*W.T
		UpdateTrm(&A22, &A21, &wrk, -1.0, 1.0, LOWER|TRANSB)

		// partially undo row pivots left of diagonal
		for k := nblk; k > 0; k-- {
			var s, d matrix.FloatMatrix
			r := p1.pivots[k-1]
			rlen := k - 1
			if r < 0 {
				r = -r
				rlen--
			}
			if r == k {
				// no pivot
				continue
			}
			ABR.SubMatrix(&s, k-1, 0, 1, rlen)
			ABR.SubMatrix(&d, r-1, 0, 1, rlen)
			Swap(&d, &s)

			if p1.pivots[k-1] < 0 {
				k-- // skip other entry in 2x2 pivots
			}
		}

		// shift pivot values
		for k, n := range p1.pivots {
			if n > 0 {
				p1.pivots[k] += ATL.Rows()
			} else {
				p1.pivots[k] -= ATL.Rows()
			}
		}

		// zero work for debuging
		W.Scale(0.0)

		// ---------------------------------------------------------

		continue3x3to2x2(
			&ATL, &ATR,
			&ABL, &ABR, &A00, &A11, &A22, A, pBOTTOMRIGHT)
		contPivot3x1to2x1(
			&pT,
			&pB, &p0, &p1, p, pBOTTOM)
	}

	// do the last part with unblocked code
	if ABR.Cols() > 0 {
		unblkDecompBKLower(&ABR, W, &pB)
		// shift pivot values
		for k, n := range pB.pivots {
			if n > 0 {
				pB.pivots[k] += ATL.Rows()
			} else {
				pB.pivots[k] -= ATL.Rows()
			}
		}
	}
	return
}

func findAndBuildBKPivotUpper(AL, AR, WL, WR *matrix.FloatMatrix, k int) (int, int) {
	var r, q int
	var rcol, qrow, src, wk, wkp1, wrow matrix.FloatMatrix

	lc := AL.Cols() - 1
	wc := WL.Cols() - 1
	lr := AL.Rows() - 1
	// Copy AR[:,lc] to WR[:,wc] and update with WL[0:]
	AL.SubMatrix(&src, 0, lc, AL.Rows(), 1)
	WL.SubMatrix(&wk, 0, wc, AL.Rows(), 1)
	src.CopyTo(&wk)
	if k > 0 {
		WR.SubMatrix(&wrow, lr, 0, 1, WR.Cols())
		//fmt.Printf("wrow: %v\n", &wrow)
		MVMult(&wk, AR, &wrow, -1.0, 1.0, NOTRANS)
		//fmt.Printf("wk after update:\n%v\n", &wk)
	}
	if AL.Rows() == 1 {
		return -1, 1
	}
	amax := math.Abs(WL.GetAt(lr, wc))

	// find max off-diagonal on first column.
	WL.SubMatrix(&rcol, 0, wc, lr, 1)
	//fmt.Printf("rcol:\n%v\n", &rcol)
	// r is row index and rmax is its absolute value
	r = IAMax(&rcol)
	rmax := math.Abs(rcol.GetAt(r, 0))
	//fmt.Printf("r=%d, amax=%e, rmax=%e\n", r, amax, rmax)
	if amax >= bkALPHA*rmax {
		// no pivoting, 1x1 diagonal
		return -1, 1
	}

	// Now we need to copy row r to WR[:,wc-1] and update it
	WL.SubMatrix(&wkp1, 0, wc-1, AL.Rows(), 1)
	if r > 0 {
		// above the diagonal part of AL
		AL.SubMatrix(&qrow, 0, r, r, 1)
		qrow.CopyTo(&wkp1)
	}
	//fmt.Printf("m(AR)=%d, r=%d, qrow: %v\n", AL.Rows(), r, &qrow)
	var wkr matrix.FloatMatrix
	AL.SubMatrix(&qrow, r, r, 1, AL.Rows()-r)
	wkp1.SubMatrix(&wkr, r, 0, AL.Rows()-r, 1)
	qrow.CopyTo(&wkr)
	//fmt.Printf("m(AR)=%d, r=%d, qrow: %v\n", AR.Rows(), r, &qrow)
	if k > 0 {
		// update wkp1
		WR.SubMatrix(&wrow, r, 0, 1, WR.Cols())
		//fmt.Printf("initial wpk1:\n%v\n", &wkp1)
		MVMult(&wkp1, AR, &wrow, -1.0, 1.0, NOTRANS)
	}
	//fmt.Printf("updated wpk1:\n%v\n", &wkp1)

	// set on-diagonal entry to zero to avoid hitting it.
	p1 := wkp1.GetAt(r, 0)
	wkp1.SetAt(r, 0, 0.0)
	// max off-diagonal on r'th column/row at index q
	q = IAMax(&wkp1)
	qmax := math.Abs(wkp1.GetAt(q, 0))
	wkp1.SetAt(r, 0, p1)
	//fmt.Printf("blk: r=%d, q=%d, amax=%e, rmax=%e, qmax=%e\n", r, q, amax, rmax, qmax)

	if amax >= bkALPHA*rmax*(rmax/qmax) {
		// no pivoting, 1x1 diagonal
		return -1, 1
	}
	// if q == r then qmax is not off-diagonal, qmax == WR[r,1] and
	// we get 1x1 pivot as following is always true
	if math.Abs(WL.GetAt(r, wc-1)) >= bkALPHA*qmax {
		// 1x1 pivoting and interchange with k, r
		// pivot row in column WR[:,1] to W[:,0]
		//p1 := WL.GetAt(r, wc-1)
		WL.SubMatrix(&src, 0, wc-1, AL.Rows(), 1)
		WL.SubMatrix(&wkp1, 0, wc, AL.Rows(), 1)
		src.CopyTo(&wkp1)
		wkp1.SetAt(-1, 0, src.GetAt(r, 0))
		wkp1.SetAt(r, 0, src.GetAt(-1, 0))
		return r, 1
	} else {
		// 2x2 pivoting and interchange with k+1, r
		return r, 2
	}
	return -1, 1
}

func unblkBoundedBKUpper(A, wrk *matrix.FloatMatrix, p *pPivots, ncol int) (error, int) {
	var err error
	var ATL, ATR, ABL, ABR matrix.FloatMatrix
	var A00, a01, A02, a11, a12t, A22, a11inv matrix.FloatMatrix
	var w00, w01, w11 matrix.FloatMatrix
	var cwrk matrix.FloatMatrix
	var wx, Ax, wz matrix.FloatMatrix
	var pT, pB, p0, p1, p2 pPivots

	err = nil
	nc := 0
	if ncol > A.Cols() {
		ncol = A.Cols()
	}

	partition2x2(
		&ATL, &ATR,
		&ABL, &ABR, A, 0, 0, pBOTTOMRIGHT)
	partitionPivot2x1(
		&pT,
		&pB, p, 0, pBOTTOM)

	// permanent working space for symmetric inverse of a11
	wrk.SubMatrix(&a11inv, wrk.Rows()-2, 0, 2, 2)
	a11inv.SetAt(0, 1, -1.0)
	a11inv.SetAt(1, 0, -1.0)

	for ATL.Cols() > 0 && nc < ncol {

		partition2x2(
			&w00, &w01,
			nil, &w11, wrk, nc, nc, pBOTTOMRIGHT)
		merge1x2(&wx, &w00, &w01)
		merge1x2(&Ax, &ATL, &ATR)

		//fmt.Printf("ATL:\n%v\n", &ATL)
		r, np := findAndBuildBKPivotUpper(&ATL, &ATR, &w00, &w01, nc)
		//fmt.Printf("[w00;w01]:\n%v\n", &wx)
		//fmt.Printf("after find: r=%d, np=%d, ncol=%d, nc=%d\n", r, np, ncol, nc)
		w00.SubMatrix(&wz, 0, w00.Cols()-2, w00.Rows(), 2)
		if np > ncol-nc {
			// next pivot does not fit into ncol columns, restore last column,
			// return with number of factorized columns
			return err, nc
		}
		if r != -1 {
			// pivoting needed; np == 1, last row; np == 2; next to last rows
			nrow := ATL.Rows() - np
			applyBKPivotSym(&ATL, nrow, r, UPPER)
			// swap left hand rows to get correct updates
			swapRows(&ATR, nrow, r)
			swapRows(&w01, nrow, r)
			if np == 2 {
				/* pivot block on diagonal; -1,-1
				 * [r, r] | [r ,-1]
				 * ----------------  2-by-2 pivot, swapping [1,0] and [r,0]
				 * [r,-1] | [-1,-1]
				 */
				t0 := w00.GetAt(-2, -1)
				tr := w00.GetAt(r, -1)
				//fmt.Printf("nc=%d, t0=%e, tr=%e\n", nc, t0, tr)
				w00.SetAt(-2, -1, tr)
				w00.SetAt(r, -1, t0)
				// interchange diagonal entries on w11[:,1]
				t0 = w00.GetAt(-2, -2)
				tr = w00.GetAt(r, -2)
				w00.SetAt(-2, -2, tr)
				w00.SetAt(r, -2, t0)
				//fmt.Printf("wrk:\n%v\n", &wz)
			}
			//fmt.Printf("pivoted A:\n%v\n", &Ax)
			//fmt.Printf("pivoted wrk:\n%v\n", &wx)
		}

		// repartition according the pivot size
		repartition2x2to3x3(&ATL,
			&A00, &a01, &A02,
			nil, &a11, &a12t,
			nil, nil, &A22 /**/, A, np, pTOPLEFT)
		repartPivot2x1to3x1(&pT,
			&p0,
			&p1,
			&p2 /**/, p, np, pTOP)
		// ------------------------------------------------------------

		wlc := w00.Cols() - np
		//wlr := w00.Rows() - 1
		w00.SubMatrix(&cwrk, 0, wlc, a01.Rows(), np)
		if np == 1 {
			//fmt.Printf("wz:\n%v\n", &wz)
			//fmt.Printf("a11 <-- %e\n", w00.GetAt(a01.Rows(), wlc))

			//w00.SubMatrix(&cwrk, 0, wlc-np+1, a01.Rows(), np)
			a11.SetAt(0, 0, w00.GetAt(a01.Rows(), wlc))
			// a21 = a21/a11
			//fmt.Printf("np == 1: pre-update a01\n%v\n", &a01)
			cwrk.CopyTo(&a01)
			InvScale(&a01, a11.Float())
			//fmt.Printf("np == 1: cwrk\n%v\na21\n%v\n", &cwrk, &a21)
			// store pivot point relative to original matrix
			if r == -1 {
				p1.pivots[0] = ATL.Rows()
			} else {
				p1.pivots[0] = r + 1
			}
		} else if np == 2 {
			/*         d | b
			 * w00 == ------
			 *         . | a
			 */
			a := w00.GetAt(-1, -1)
			b := w00.GetAt(-2, -1)
			d := w00.GetAt(-2, -2)
			a11inv.SetAt(1, 1, d/b)
			a11inv.SetAt(0, 0, a/b)
			// denominator: (a/b)*(d/b)-1.0 == (a*d - b^2)/b^2
			scale := 1.0 / ((a/b)*(d/b) - 1.0)
			scale /= b
			//fmt.Printf("a11inv:\n%v\n", &a11inv)

			// a01 = a01*a11.-1
			Mult(&a01, &cwrk, &a11inv, scale, 0.0, NOTRANS)
			a11.SetAt(1, 1, a)
			a11.SetAt(0, 1, b)
			a11.SetAt(0, 0, d)

			// store pivot point relative to original matrix
			p1.pivots[0] = -(r + 1)
			p1.pivots[1] = p1.pivots[0]
		}

		//fmt.Printf("end-of-loop: Ax r=%d, nc=%d. np=%d\n%v\n", r, nc, np, &Ax)
		//fmt.Printf("wx m(wblk)=%d:\n%v\n", m(&wx), &wx)

		// ------------------------------------------------------------
		nc += np
		continue3x3to2x2(
			&ATL, &ATR,
			&ABL, &ABR, &A00, &a11, &A22, A, pTOPLEFT)
		contPivot3x1to2x1(
			&pT,
			&pB, &p0, &p1, p, pTOP)

	}
	return err, nc
}

func blkDecompBKUpper(A, W *matrix.FloatMatrix, p *pPivots, nb int) (err error) {
	var ATL, ATR, ABL, ABR matrix.FloatMatrix
	var A00, A01, A02, A11, A12, A22 matrix.FloatMatrix
	var wrk matrix.FloatMatrix
	var pT, pB, p0, p1, p2 pPivots
	var nblk int = 0

	err = nil
	partition2x2(
		&ATL, &ATR,
		&ABL, &ABR, A, 0, 0, pBOTTOMRIGHT)
	partitionPivot2x1(
		&pT,
		&pB, p, 0, pBOTTOM)

	for ATL.Cols() >= nb {
		err, nblk = unblkBoundedBKUpper(&ATL, W, &pT, nb)

		// repartition nblk size
		repartition2x2to3x3(&ATL,
			&A00, &A01, &A02,
			nil, &A11, &A12,
			nil, nil, &A22, A, nblk, pTOPLEFT)
		repartPivot2x1to3x1(&pT,
			&p0, &p1, &p2 /**/, p, nblk, pTOP)

		// --------------------------------------------------------
		// here [A01;A11] has been decomposed by unblkBoundedBKUpper()
		// Now we need update A00

		// wrk is original A01; D1*L01.T
		W.SubMatrix(&wrk, 0, W.Cols()-nblk, A01.Rows(), nblk)

		// A00 = A00 - L01*D1*L01.T = A00 - L01*W.T
		UpdateTrm(&A00, &A01, &wrk, -1.0, 1.0, UPPER|TRANSB)

		// partially undo row pivots right of diagonal
		for k := 0; k < nblk; k++ {
			var s, d matrix.FloatMatrix
			r := p1.pivots[k]
			colno := A00.Cols() + k
			np := 1
			if r < 0 {
				r = -r
				np = 2
			}
			rlen := ATL.Cols() - colno - np
			//fmt.Printf("undo: k=%d, r=%d, colno=%d, rlen=%d\n", k, r, colno, rlen)
			if r == colno+1 {
				// no pivot
				continue
			}
			ATL.SubMatrix(&s, colno, colno+np, 1, rlen)
			ATL.SubMatrix(&d, r-1, colno+np, 1, rlen)
			//fmt.Printf("s %d: %v\n", colno, &s)
			//fmt.Printf("d %d: %v\n", r-1,   &d)
			Swap(&d, &s)

			if p1.pivots[k] < 0 {
				k++ // skip other entry in 2x2 pivots
			}
		}

		// ---------------------------------------------------------

		continue3x3to2x2(
			&ATL, &ATR,
			&ABL, &ABR, &A00, &A11, &A22, A, pTOPLEFT)
		contPivot3x1to2x1(
			&pT,
			&pB, &p0, &p1, p, pTOP)
	}

	// do the last part with unblocked code
	if ATL.Cols() > 0 {
		unblkDecompBKUpper(&ATL, W, &pT)
	}
	return
}

func DecomposeBK(A, W *matrix.FloatMatrix, ipiv []int, flags Flags, nb int) (*matrix.FloatMatrix, error) {
	var err error = nil
	for k, _ := range ipiv {
		ipiv[k] = 0
	}
	if A.Cols() < nb || nb == 0 {
		if W.Cols() < 2 || W.Rows() < A.Rows() {
			return nil, errors.New("Workspace too small")
		}
		if flags&LOWER != 0 {
			err, _ = unblkDecompBKLower(A, W, &pPivots{ipiv})
		} else if flags&UPPER != 0 {
			err, _ = unblkDecompBKUpper(A, W, &pPivots{ipiv})
		}
	} else {
		if W.Cols() < nb+1 || W.Rows() < A.Rows() {
			return nil, errors.New("Workspace too small")
		}
		if flags&LOWER != 0 {
			err = blkDecompBKLower(A, W, &pPivots{ipiv}, nb)
		} else if flags&UPPER != 0 {
			err = blkDecompBKUpper(A, W, &pPivots{ipiv}, nb)
		}
	}
	return A, err
}

// Local Variables:
// tab-width: 4
// indent-tabs-mode: nil
// End:
