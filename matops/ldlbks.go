// Copyright (c) Harri Rautila, 2012,2013

// This file is part of github.com/hrautila/matops package. It is free software,
// distributed under the terms of GNU Lesser General Public License Version 3, or
// any later version. See the COPYING tile included in this archive.

package matops

import (
	"fmt"
	"github.com/henrylee2cn/algorithm/matrix"
)

func unblkSolveBKLower(B, A /*, wrk*/ *matrix.FloatMatrix, p *pPivots, phase int) error {
	var err error
	var ATL, ATR, ABL, ABR matrix.FloatMatrix
	var A00, a10t, a11, A20, a21, A22 /*, a11inv*/ matrix.FloatMatrix
	var Aref *matrix.FloatMatrix
	var BT, BB, B0, b1, B2, Bx matrix.FloatMatrix
	var pT, pB, p0, p1, p2 pPivots
	var aStart, aDir, bStart, bDir pDirection
	var nc int

	err = nil
	np := 0

	if phase == 1 {
		aStart = pTOPLEFT
		aDir = pBOTTOMRIGHT
		bStart = pTOP
		bDir = pBOTTOM
		nc = 1
		Aref = &ABR
	} else {
		aStart = pBOTTOMRIGHT
		aDir = pTOPLEFT
		bStart = pBOTTOM
		bDir = pTOP
		nc = A.Rows()
		Aref = &ATL
	}
	partition2x2(
		&ATL, &ATR,
		&ABL, &ABR, A, 0, 0, aStart)
	partition2x1(
		&BT,
		&BB, B, 0, bStart)
	partitionPivot2x1(
		&pT,
		&pB, p, 0, bStart)

	// ABR.Cols() == 0 is end of matrix,
	for Aref.Cols() > 0 {

		// see if next diagonal block is 1x1 or 2x2
		np = 1
		if p.pivots[nc-1] < 0 {
			np = 2
		}
		//fmt.Printf("nc=%d, np=%d, m(ABR)=%d\n", nc, np, m(&ABR))

		// repartition according the pivot size
		repartition2x2to3x3(&ATL,
			&A00, nil, nil,
			&a10t, &a11, nil,
			&A20, &a21, &A22 /**/, A, np, aDir)
		repartition2x1to3x1(&BT,
			&B0,
			&b1,
			&B2 /**/, B, np, bDir)
		repartPivot2x1to3x1(&pT,
			&p0,
			&p1,
			&p2 /**/, p, np, bDir)
		// ------------------------------------------------------------

		switch phase {
		case 1:
			// computes D.-1*(L.-1*B)
			if np == 1 {
				if p1.pivots[0] != nc {
					// swap rows on bottom part of B
					//fmt.Printf("1x1 pivot top with %d [%d]\n", p1.pivots[0], p1.pivots[0]-BT.Rows())
					swapRows(&BB, 0, p1.pivots[0]-BT.Rows()-1)
				}
				// B2 = B2 - a21*b1
				MVRankUpdate(&B2, &a21, &b1, -1.0)
				// b1 = b1/d1
				InvScale(&b1, a11.Float())
				nc += 1
			} else if np == 2 {
				if p1.pivots[0] != -nc {
					// swap rows on bottom part of B
					//fmt.Printf("2x2 pivot %d with %d [%d]\n", nc+1, -p1.pivots[0])
					//fmt.Printf("pre :\n%v\n", B)
					swapRows(&BB, 1, -p1.pivots[0]-BT.Rows()-1)
					//fmt.Printf("post:\n%v\n", B)
				}
				b := a11.GetAt(1, 0)
				apb := a11.GetAt(0, 0) / b
				dpb := a11.GetAt(1, 1) / b
				// (a/b)*(d/b)-1.0 == (a*d - b^2)/b^2
				scale := apb*dpb - 1.0
				scale *= b

				// B2 = B2 - a21*b1
				Mult(&B2, &a21, &b1, -1.0, 1.0, NOTRANS)
				// b1 = a11.-1*b1.T
				//(2x2 block, no subroutine for doing this in-place)
				for k := 0; k < b1.Cols(); k++ {
					s0 := b1.GetAt(0, k)
					s1 := b1.GetAt(1, k)
					b1.SetAt(0, k, (dpb*s0-s1)/scale)
					b1.SetAt(1, k, (apb*s1-s0)/scale)
				}
				nc += 2
			}
		case 2:
			if np == 1 {
				MVMult(&b1, &B2, &a21, -1.0, 1.0, TRANSA)
				if p1.pivots[0] != nc {
					// swap rows on bottom part of B
					//fmt.Printf("1x1 pivot top with %d [%d]\n", p1.pivots[0], p1.pivots[0]-BT.Rows())
					merge2x1(&Bx, &b1, &B2)
					swapRows(&Bx, 0, p1.pivots[0]-BT.Rows())
				}
				nc -= 1
			} else if np == 2 {
				Mult(&b1, &a21, &B2, -1.0, 1.0, TRANSA)
				if p1.pivots[0] != -nc {
					// swap rows on bottom part of B
					//fmt.Printf("2x2 pivot %d with %d\n", nc, -p1.pivots[0])
					merge2x1(&Bx, &b1, &B2)
					//fmt.Printf("pre :\n%v\n", B)
					swapRows(&Bx, 1, -p1.pivots[0]-BT.Rows()+1)
					//fmt.Printf("post:\n%v\n", B)
				}
				nc -= 2
			}
		}

		// ------------------------------------------------------------

		continue3x3to2x2(
			&ATL, &ATR,
			&ABL, &ABR, &A00, &a11, &A22, A, aDir)
		continue3x1to2x1(
			&BT,
			&BB, &B0, &b1, B, bDir)
		contPivot3x1to2x1(
			&pT,
			&pB, &p0, &p1, p, bDir)

	}
	return err
}

func unblkSolveBKUpper(B, A *matrix.FloatMatrix, p *pPivots, phase int) error {
	var err error
	var ATL, ATR, ABL, ABR matrix.FloatMatrix
	var A00, a01, A02, a11, a12t, A22 matrix.FloatMatrix
	var Aref *matrix.FloatMatrix
	var BT, BB, B0, b1, B2, Bx matrix.FloatMatrix
	var pT, pB, p0, p1, p2 pPivots
	var aStart, aDir, bStart, bDir pDirection
	var nc int

	err = nil
	np := 0

	if phase == 2 {
		aStart = pTOPLEFT
		aDir = pBOTTOMRIGHT
		bStart = pTOP
		bDir = pBOTTOM
		nc = 1
		Aref = &ABR
	} else {
		aStart = pBOTTOMRIGHT
		aDir = pTOPLEFT
		bStart = pBOTTOM
		bDir = pTOP
		nc = A.Rows()
		Aref = &ATL
	}
	partition2x2(
		&ATL, &ATR,
		&ABL, &ABR, A, 0, 0, aStart)
	partition2x1(
		&BT,
		&BB, B, 0, bStart)
	partitionPivot2x1(
		&pT,
		&pB, p, 0, bStart)

	// ABR.Cols() == 0 is end of matrix,
	for Aref.Cols() > 0 {

		// see if next diagonal block is 1x1 or 2x2
		np = 1
		if p.pivots[nc-1] < 0 {
			np = 2
		}
		fmt.Printf("nc=%d, np=%d, m(ABR)=%d\n", nc, np, m(&ABR))

		// repartition according the pivot size
		repartition2x2to3x3(&ATL,
			&A00, &a01, &A02,
			nil, &a11, &a12t,
			nil, nil, &A22 /**/, A, np, aDir)
		repartition2x1to3x1(&BT,
			&B0,
			&b1,
			&B2 /**/, B, np, bDir)
		repartPivot2x1to3x1(&pT,
			&p0,
			&p1,
			&p2 /**/, p, np, bDir)
		// ------------------------------------------------------------

		switch phase {
		case 1:
			// computes D.-1*(L.-1*B)
			if np == 1 {
				if p1.pivots[0] != nc {
					// swap rows in top part of B
					//fmt.Printf("1x1 pivot top with %d [%d]\n", p1.pivots[0], p1.pivots[0]-BT.Rows())
					swapRows(&BT, BT.Rows()-1, p1.pivots[0]-1)
				}
				// B2 = B2 - a21*b1
				MVRankUpdate(&B2, &a01, &b1, -1.0)
				// b1 = b1/d1
				InvScale(&b1, a11.Float())
				nc += 1
			} else if np == 2 {
				if p1.pivots[0] != -nc {
					// swap rows on bottom part of B
					//fmt.Printf("2x2 pivot %d with %d [%d]\n", nc+1, -p1.pivots[0])
					//fmt.Printf("pre :\n%v\n", B)
					swapRows(&BT, BT.Rows()-2, -p1.pivots[0]-1)
					//fmt.Printf("post:\n%v\n", B)
				}
				b := a11.GetAt(0, 1)
				apb := a11.GetAt(0, 0) / b
				dpb := a11.GetAt(1, 1) / b
				// (a/b)*(d/b)-1.0 == (a*d - b^2)/b^2
				scale := apb*dpb - 1.0
				scale *= b

				// B2 = B2 - a21*b1
				Mult(&B2, &a01, &b1, -1.0, 1.0, NOTRANS)
				// b1 = a11.-1*b1.T
				//(2x2 block, no subroutine for doing this in-place)
				for k := 0; k < b1.Cols(); k++ {
					s0 := b1.GetAt(0, k)
					s1 := b1.GetAt(1, k)
					b1.SetAt(0, k, (dpb*s0-s1)/scale)
					b1.SetAt(1, k, (apb*s1-s0)/scale)
				}
				nc += 2
			}
		case 2:
			if np == 1 {
				MVMult(&b1, &B2, &a01, -1.0, 1.0, TRANSA)
				if p1.pivots[0] != nc {
					// swap rows on bottom part of B
					//fmt.Printf("1x1 pivot top with %d [%d]\n", p1.pivots[0], p1.pivots[0]-BT.Rows())
					merge2x1(&Bx, &B0, &b1)
					swapRows(&Bx, Bx.Rows()-1, p1.pivots[0]-1)
				}
				nc -= 1
			} else if np == 2 {
				Mult(&b1, &a01, &B2, -1.0, 1.0, TRANSA)
				if p1.pivots[0] != -nc {
					// swap rows on bottom part of B
					//fmt.Printf("2x2 pivot %d with %d\n", nc, -p1.pivots[0])
					merge2x1(&Bx, &B0, &b1)
					//fmt.Printf("pre :\n%v\n", B)
					swapRows(&Bx, Bx.Rows()-2, -p1.pivots[0]-1)
					//fmt.Printf("post:\n%v\n", B)
				}
				nc -= 2
			}
		}

		// ------------------------------------------------------------

		continue3x3to2x2(
			&ATL, &ATR,
			&ABL, &ABR, &A00, &a11, &A22, A, aDir)
		continue3x1to2x1(
			&BT,
			&BB, &B0, &b1, B, bDir)
		contPivot3x1to2x1(
			&pT,
			&pB, &p0, &p1, p, bDir)

	}
	return err
}

func SolveBK(B, A *matrix.FloatMatrix, ipiv []int, flags Flags) {
	if flags&LOWER != 0 {
		// first part: Z = D.-1*(L.-1*B)
		unblkSolveBKLower(B, A, &pPivots{ipiv}, 1)
		// second part: X = L.-T*Z
		unblkSolveBKLower(B, A, &pPivots{ipiv}, 2)
	} else if flags&UPPER != 0 {
		// first part: Z = D.-1*(U.-1*B)
		unblkSolveBKUpper(B, A, &pPivots{ipiv}, 1)
		// second part: X = U.-T*Z
		unblkSolveBKUpper(B, A, &pPivots{ipiv}, 2)
	}
}

// Local Variables:
// tab-width: 4
// indent-tabs-mode: nil
// End:
