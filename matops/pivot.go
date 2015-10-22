// Copyright (c) Harri Rautila, 2012,2013

// This file is part of github.com/hrautila/matops package. It is free software,
// distributed under the terms of GNU Lesser General Public License Version 3, or
// any later version. See the COPYING tile included in this archive.

package matops

import (
	"github.com/henrylee2cn/algorithm/matrix"
	"math"
	//"fmt"
)

type pPivots struct {
	pivots []int
}

const (
	FORWARD  = 1
	BACKWARD = 2
)

func swapRows(A *matrix.FloatMatrix, src, dst int) {
	var r0, r1 matrix.FloatMatrix
	if src == dst || A.Rows() == 0 {
		return
	}
	A.SubMatrix(&r0, src, 0, 1, A.Cols())
	A.SubMatrix(&r1, dst, 0, 1, A.Cols())
	Swap(&r0, &r1)
}

func swapCols(A *matrix.FloatMatrix, src, dst int) {
	var c0, c1 matrix.FloatMatrix
	if src == dst || A.Rows() == 0 {
		return
	}
	A.SubMatrix(&c0, 0, src, A.Rows(), 1)
	A.SubMatrix(&c1, 0, dst, A.Rows(), 1)
	Swap(&c0, &c1)
}

func scalePivots(p *pPivots, offset int) {
	for k := 0; k < len(p.pivots); k++ {
		if p.pivots[k] > 0 {
			p.pivots[k] += offset
		}
	}
}

func applyPivots(A *matrix.FloatMatrix, p *pPivots) {
	for k, n := range p.pivots {
		if n > 0 {
			swapRows(A, n, k)
		}
	}
}

func applyRowPivots(A *matrix.FloatMatrix, p *pPivots, offset, dir int) {
	if dir == FORWARD {
		for k, n := range p.pivots {
			if n > 0 {
				swapRows(A, n-1-offset, k)
			}
		}
	} else if dir == BACKWARD {
		//
		for k := len(p.pivots) - 1; k >= 0; k-- {
			if p.pivots[k] > 0 {
				swapRows(A, p.pivots[k]-1-offset, k)
			}
		}
	}
}

func applyColPivots(A *matrix.FloatMatrix, p *pPivots, offset, dir int) {
	if dir == FORWARD {
		for k, n := range p.pivots {
			if n > 0 {
				swapCols(A, n-1-offset, k)
			}
		}
	} else if dir == BACKWARD {
		//
		for k := len(p.pivots) - 1; k >= 0; k-- {
			if p.pivots[k] > 0 {
				swapCols(A, p.pivots[k]-1-offset, k)
			}
		}
	}
}

// Find largest absolute value on column
func pivotIndex(A *matrix.FloatMatrix, p *pPivots) {
	max := math.Abs(A.GetAt(0, 0))
	for k := 1; k < A.Rows(); k++ {
		v := math.Abs(A.GetAt(k, 0))
		if v > max {
			p.pivots[0] = k
			max = v
		}
	}
}

/*
 * Apply diagonal pivot (row and column swapped) to symmetric matrix blocks.
 * AR[0,0] is on diagonal and AL is block to the left of diagonal and AR the
 * triangular diagonal block. Need to swap row and column.
 *
 * LOWER triangular; moving from top-left to bottom-right
 *
 *    d
 *    x  d
 *    x  x  d  |
 *    --------------------------
 *    S1 S1 S1 | P1 x  x  x  P2     -- current row
 *    x  x  x  | S2 d  x  x  x
 *    x  x  x  | S2 x  d  x  x
 *    x  x  x  | S2 x  x  d  x
 *    D1 D1 D1 | P2 D2 D2 D2 P3     -- swap with row 'index'
 *    x  x  x  | S3 x  x  x  D3 d
 *    x  x  x  | S3 x  x  x  D3 x d
 *       (ABL)          (ABR)
 *
 * UPPER triangular; moving from bottom-right to top-left
 *
 *         (ATL)             (ATR)
 *    d  x  x  D3 x  x  x | S3 x  x
 *       d  x  D3 x  x  x | S3 x  x
 *          d  D3 x  x  x | S3 x  x
 *             P3 D2 D2 D2| P2 D1 D1
 *                d  x  x | S2 x  x
 *                   d  x | S2 x  x
 *                      d | S2 x  x
 *    -----------------------------
 *                        | P1 S1 S1
 *                        |    d  x
 *                        |       d
 *                           (ABR)
 */
func applyPivotSym(AL, AR *matrix.FloatMatrix, index int, flags Flags) {
	var s, d matrix.FloatMatrix
	if flags&LOWER != 0 {
		// AL is [ABL]; AR is [ABR]; P1 is AR[0,0], P2 is AR[index, 0]
		// S1 -- D1
		AL.SubMatrix(&s, 0, 0, 1, AL.Cols())
		AL.SubMatrix(&d, index, 0, 1, AL.Cols())
		Swap(&s, &d)
		// S2 -- D2
		AR.SubMatrix(&s, 1, 0, index-1, 1)
		AR.SubMatrix(&d, index, 1, 1, index-1)
		Swap(&s, &d)
		// S3 -- D3
		AR.SubMatrix(&s, index+1, 0, AR.Rows()-index-1, 1)
		AR.SubMatrix(&d, index+1, index, AR.Rows()-index-1, 1)
		Swap(&s, &d)
		// swap P1 and P3
		p1 := AR.GetAt(0, 0)
		p3 := AR.GetAt(index, index)
		AR.SetAt(0, 0, p3)
		AR.SetAt(index, index, p1)
		return
	}
	if flags&UPPER != 0 {
		// AL is merged from [ATL, ATR], AR is [ABR]; P1 is AR[0, 0]; P2 is AL[index, -1]
		colno := AL.Cols() - AR.Cols()
		// S1 -- D1; S1 is on the first row of AR
		AR.SubMatrix(&s, 0, 1, 1, AR.Cols()-1)
		AL.SubMatrix(&d, index, colno+1, 1, s.Cols())
		Swap(&s, &d)
		// S2 -- D2
		AL.SubMatrix(&s, index+1, colno, AL.Rows()-index-2, 1)
		AL.SubMatrix(&d, index, index+1, 1, colno-index-1)
		Swap(&s, &d)
		// S3 -- D3
		AL.SubMatrix(&s, 0, index, index, 1)
		AL.SubMatrix(&d, 0, colno, index, 1)
		Swap(&s, &d)
		//fmt.Printf("3, AR=%v\n", AR)
		// swap P1 and P3
		p1 := AR.GetAt(0, 0)
		p3 := AL.GetAt(index, index)
		AR.SetAt(0, 0, p3)
		AL.SetAt(index, index, p1)
		return
	}
}

/*
 * Apply diagonal pivot (row and column swapped) to symmetric matrix blocks.
 * AR[0,0] is on diagonal and AL is block to the left of diagonal and AR the
 * triangular diagonal block. Need to swap row and column.
 *
 * LOWER triangular; moving from top-left to bottom-right
 *
 *    d
 *    x  d |
 *    --------------------------
 *    x  x | d
 *    S1 S1| S1 P1 x  x  x  P2     -- current row/col 'srcix'
 *    x  x | x  S2 d  x  x  x
 *    x  x | x  S2 x  d  x  x
 *    x  x | x  S2 x  x  d  x
 *    D1 D1| D1 P2 D2 D2 D2 P3     -- swap with row/col 'dstix'
 *    x  x | x  S3 x  x  x  D3 d
 *    x  x | x  S3 x  x  x  D3 x d
 *    (ABL)          (ABR)
 *
 * UPPER triangular; moving from bottom-right to top-left
 *
 *         (ATL)                  (ATR)
 *    d  x  x  D3 x  x  x  S3 x | x
 *       d  x  D3 x  x  x  S3 x | x
 *          d  D3 x  x  x  S3 x | x
 *             P3 D2 D2 D2 P2 D1| D1  -- dstinx
 *                d  x  x  S2 x | x
 *                   d  x  S2 x | x
 *                      d  S2 x | x
 *                         P1 S1| S1  -- srcinx
 *                            d | x
 *    -----------------------------
 *                              | d
 *                           (ABR)
 */
func applyPivotSym2(AL, AR *matrix.FloatMatrix, srcix, dstix int, flags Flags) {
	var s, d matrix.FloatMatrix
	if flags&LOWER != 0 {
		// AL is [ABL]; AR is [ABR]; P1 is AR[0,0], P2 is AR[index, 0]
		// S1 -- D1
		AL.SubMatrix(&s, srcix, 0, 1, AL.Cols())
		AL.SubMatrix(&d, dstix, 0, 1, AL.Cols())
		Swap(&s, &d)
		if srcix > 0 {
			AR.SubMatrix(&s, srcix, 0, 1, srcix)
			AR.SubMatrix(&d, dstix, 0, 1, srcix)
			Swap(&s, &d)
		}
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
		// S1 -- D1;
		AR.SubMatrix(&s, srcix, 0, 1, AR.Cols())
		AR.SubMatrix(&d, dstix, 0, 1, AR.Cols())
		Swap(&s, &d)
		if srcix < AL.Cols()-1 {
			// not the corner element
			AL.SubMatrix(&s, srcix, srcix+1, 1, srcix)
			AL.SubMatrix(&d, dstix, srcix+1, 1, srcix)
			Swap(&s, &d)
		}
		// S2 -- D2
		AL.SubMatrix(&s, dstix+1, srcix, srcix-dstix-1, 1)
		AL.SubMatrix(&d, dstix, dstix+1, 1, srcix-dstix-1)
		Swap(&s, &d)
		// S3 -- D3
		AL.SubMatrix(&s, 0, srcix, dstix, 1)
		AL.SubMatrix(&d, 0, dstix, dstix, 1)
		Swap(&s, &d)
		//fmt.Printf("3, AR=%v\n", AR)
		// swap P1 and P3
		p1 := AR.GetAt(0, 0)
		p3 := AL.GetAt(dstix, dstix)
		AR.SetAt(srcix, srcix, p3)
		AL.SetAt(dstix, dstix, p1)
		return
	}
}

func ApplyRowPivots(A *matrix.FloatMatrix, ipiv []int, direction int) {
	p := &pPivots{ipiv}
	applyRowPivots(A, p, 0, direction)
}

func NumPivots(ipiv []int) int {
	count := 0
	for _, n := range ipiv {
		if n != 0 {
			count += 1
		}
	}
	return count
}

// Local Variables:
// tab-width: 4
// indent-tabs-mode: nil
// End:
