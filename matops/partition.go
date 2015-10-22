// Copyright (c) Harri Rautila, 2012,2013

// This file is part of github.com/hrautila/matops package. It is free software,
// distributed under the terms of GNU Lesser General Public License Version 3, or
// any later version. See the COPYING tile included in this archive.

package matops

import (
	"github.com/henrylee2cn/algorithm/matrix"
	//"fmt"
)

// Functions here are support functions for libFLAME-like implementation
// of various linear algebra algorithms.

type pDirection int

const (
	pLEFT = iota
	pRIGHT
	pTOP
	pBOTTOM
	pTOPLEFT
	pBOTTOMRIGHT
)

/*
 Partition p to 2 by 1 blocks.

        AT
  A --> --
        AB

 Parameter nb is initial block size for AT (pTOP) or AB (pBOTTOM).
*/
func partition2x1(AT, AB, A *matrix.FloatMatrix, nb int, side pDirection) {
	if nb > A.Rows() {
		nb = A.Rows()
	}
	switch side {
	case pTOP:
		A.SubMatrix(AT, 0, 0, nb, A.Cols())
		A.SubMatrix(AB, nb, 0, A.Rows()-nb, A.Cols())
	case pBOTTOM:
		A.SubMatrix(AT, 0, 0, A.Rows()-nb, A.Cols())
		A.SubMatrix(AB, A.Rows()-nb, 0, nb, A.Cols())
	}
}

/*
 Repartition 2 by 1 block to 3 by 1 block.

           AT      A0            AT       A0
 pBOTTOM: --  --> --   ; pTOP:   --  -->  A1
           AB      A1            AB       --
                   A2                     A2

*/
func repartition2x1to3x1(AT, A0, A1, A2, A *matrix.FloatMatrix, nb int, pdir pDirection) {
	nT := AT.Rows()
	switch pdir {
	case pBOTTOM:
		if nT+nb > A.Rows() {
			nb = A.Rows() - nT
		}
		A.SubMatrix(A0, 0, 0, nT, A.Cols())
		A.SubMatrix(A1, nT, 0, nb, A.Cols())
		A.SubMatrix(A2, nT+nb, 0, A.Rows()-nT-nb, A.Cols())
	case pTOP:
		if nT < nb {
			nb = nT
		}
		A.SubMatrix(A0, 0, 0, nT-nb, A.Cols())
		A.SubMatrix(A1, nT-nb, 0, nb, A.Cols())
		A.SubMatrix(A2, nT, 0, A.Rows()-nT, A.Cols())
	}
}

/*
 Continue with 2 by 1 block from 3 by 1 block.

           AT      A0            AT       A0
 pBOTTOM: --  <--  A1   ; pTOP:   -- <--  --
           AB      --            AB       A1
                   A2                     A2

*/
func continue3x1to2x1(AT, AB, A0, A1, A *matrix.FloatMatrix, pdir pDirection) {
	n0 := A0.Rows()
	n1 := A1.Rows()
	switch pdir {
	case pBOTTOM:
		A.SubMatrix(AT, 0, 0, n0+n1, A.Cols())
		A.SubMatrix(AB, n0+n1, 0, A.Rows()-n0-n1, A.Cols())
	case pTOP:
		A.SubMatrix(AT, 0, 0, n0, A.Cols())
		A.SubMatrix(AB, n0, 0, A.Rows()-n0, A.Cols())
	}
}

/*
 * Merge 1 by 1 block from 2 by 1 block.
 *
 *          AT
 * Abkl <-- --
 *          AB
 *
 */
func merge2x1(ABLK, AT, AB *matrix.FloatMatrix) {
	AT.SubMatrix(ABLK, 0, 0, AT.Rows()+AB.Rows(), AT.Cols())
}

/*
 * Merge 1 by 1 block from 1 by 2 block.
 *
 * ABLK <--  AL | AR
 *
 */
func merge1x2(ABLK, AL, AR *matrix.FloatMatrix) {
	AL.SubMatrix(ABLK, 0, 0, AL.Rows(), AL.Cols()+AR.Cols())
}

/*
 Partition A to 1 by 2 blocks.

  A -->  AL | AR

 Parameter nb is initial block size for AL (pLEFT) or AR (pRIGHT).
*/
func partition1x2(AL, AR, A *matrix.FloatMatrix, nb int, side pDirection) {
	if nb > A.Cols() {
		nb = A.Cols()
	}
	switch side {
	case pLEFT:
		A.SubMatrix(AL, 0, 0, A.Rows(), nb)
		A.SubMatrix(AR, 0, nb, A.Rows(), A.Cols()-nb)
	case pRIGHT:
		A.SubMatrix(AL, 0, 0, A.Rows(), A.Cols()-nb)
		A.SubMatrix(AR, 0, A.Cols()-nb, A.Rows(), nb)
	}
}

/*
 Repartition 1 by 2 blocks to 1 by 3 blocks.

 pRIGHT: AL | AR  -->  A0 | A1 A2
 pLEFT:  AL | AR  -->  A0 A1 | A2

 Parameter As is left or right block of original 1x2 block.
*/
func repartition1x2to1x3(AL, A0, A1, A2, A *matrix.FloatMatrix, nb int, pdir pDirection) {
	k := AL.Cols()
	switch pdir {
	case pRIGHT:
		if k+nb > A.Cols() {
			nb = A.Cols() - k
		}
		// A0 is AL; [A1; A2] is AR
		A.SubMatrix(A0, 0, 0, A.Rows(), k)
		A.SubMatrix(A1, 0, k, A.Rows(), nb)
		A.SubMatrix(A2, 0, k+nb, A.Rows(), A.Cols()-nb-k)
	case pLEFT:
		if nb > k {
			nb = k
		}
		// A2 is AR; [A0; A1] is AL
		A.SubMatrix(A0, 0, 0, A.Rows(), k-nb)
		A.SubMatrix(A1, 0, k-nb, A.Rows(), nb)
		A.SubMatrix(A2, 0, k, A.Rows(), A.Cols()-k)
	}
}

/*
 Repartition 1 by 2 blocks to 1 by 3 blocks.

 pRIGHT: AL | AR  --  A0 A1 | A2
 pLEFT:  AL | AR  <--  A0 | A1 A2

*/
func continue1x3to1x2(AL, AR, A0, A1, A *matrix.FloatMatrix, pdir pDirection) {

	k := A0.Cols()
	nb := A1.Cols()
	switch pdir {
	case pRIGHT:
		// AL is [A0; A1], AR is A2
		A.SubMatrix(AL, 0, 0, A.Rows(), k+nb)
		A.SubMatrix(AR, 0, AL.Cols(), A.Rows(), A.Cols()-AL.Cols())
	case pLEFT:
		// AL is A0; AR is [A1; A2]
		A.SubMatrix(AL, 0, 0, A.Rows(), k)
		A.SubMatrix(AR, 0, k, A.Rows(), A.Cols()-k)
	}
}

/*
 Partition A to 2 by 2 blocks.

           ATL | ATR
  A  -->   =========
           ABL | ABR

 Parameter nb is initial block size for ATL in column direction and mb in row direction.
 ATR and ABL may be nil pointers.
*/
func partition2x2(ATL, ATR, ABL, ABR, A *matrix.FloatMatrix, mb, nb int, side pDirection) {
	switch side {
	case pTOPLEFT:
		A.SubMatrix(ATL, 0, 0, mb, nb)
		if ATR != nil {
			A.SubMatrix(ATR, 0, nb, mb, A.Cols()-nb)
		}
		if ABL != nil {
			A.SubMatrix(ABL, mb, 0, A.Rows()-mb, nb)
		}
		A.SubMatrix(ABR, mb, nb)
	case pBOTTOMRIGHT:
		A.SubMatrix(ATL, 0, 0, A.Rows()-mb, A.Cols()-nb)
		if ATR != nil {
			A.SubMatrix(ATR, 0, A.Cols()-nb, A.Rows()-mb, nb)
		}
		if ABL != nil {
			A.SubMatrix(ABL, A.Rows()-mb, 0, mb, nb)
		}
		A.SubMatrix(ABR, A.Rows()-mb, A.Cols()-nb)
	}
}

/*
 Repartition 2 by 2 blocks to 3 by 3 blocks.

                      A00 | A01 : A02
   ATL | ATR   nb     ===============
   =========   -->    A10 | A11 : A12
   ABL | ABR          ---------------
                      A20 | A21 : A22

 ATR, ABL, ABR implicitely defined by ATL and A.
 It is valid to have either the strictly upper or lower submatrices as nil values.

*/
func repartition2x2to3x3(ATL,
	A00, A01, A02, A10, A11, A12, A20, A21, A22, A *matrix.FloatMatrix, nb int, pdir pDirection) {

	k := ATL.Rows()
	switch pdir {
	case pBOTTOMRIGHT:
		if k+nb > A.Cols() {
			nb = A.Cols() - k
		}
		A.SubMatrix(A00, 0, 0, k, k)
		if A01 != nil {
			A.SubMatrix(A01, 0, k, k, nb)
		}
		if A02 != nil {
			A.SubMatrix(A02, 0, k+nb, k, A.Cols()-k-nb)
		}

		if A10 != nil {
			A.SubMatrix(A10, k, 0, nb, k)
		}
		A.SubMatrix(A11, k, k, nb, nb)
		if A12 != nil {
			A.SubMatrix(A12, k, k+nb, nb, A.Cols()-k-nb)
		}

		if A20 != nil {
			A.SubMatrix(A20, k+nb, 0, A.Rows()-k-nb, k)
		}
		if A21 != nil {
			A.SubMatrix(A21, k+nb, k, A.Rows()-k-nb, nb)
		}
		A.SubMatrix(A22, k+nb, k+nb)
	case pTOPLEFT:
		if nb > k {
			nb = k
		}
		// move towards top left corner
		A.SubMatrix(A00, 0, 0, k-nb, k-nb)
		if A01 != nil {
			A.SubMatrix(A01, 0, k-nb, k-nb, nb)
		}
		if A02 != nil {
			A.SubMatrix(A02, 0, k, k-nb, A.Cols()-k)
		}

		if A10 != nil {
			A.SubMatrix(A10, k-nb, 0, nb, k-nb)
		}
		A.SubMatrix(A11, k-nb, k-nb, nb, nb)
		if A12 != nil {
			A.SubMatrix(A12, k-nb, k, nb, A.Cols()-k)
		}

		if A20 != nil {
			A.SubMatrix(A20, k, 0, A.Rows()-k, k-nb)
		}
		if A21 != nil {
			A.SubMatrix(A21, k, k-nb, A.Rows()-k, nb)
		}
		A.SubMatrix(A22, k, k)
	}
}

/*
 Redefine 2 by 2 blocks from 3 by 3 partition.

                      A00 : A01 | A02
   ATL | ATR   nb     ---------------
   =========   <--    A10 : A11 | A12
   ABL | ABR          ===============
                      A20 : A21 | A22

 New division of ATL, ATR, ABL, ABR defined by diagonal entries A00, A11, A22
*/
func continue3x3to2x2(
	ATL, ATR, ABL, ABR,
	A00, A11, A22, A *matrix.FloatMatrix, pdir pDirection) {

	k := A00.Rows()
	mb := A11.Cols()
	switch pdir {
	case pBOTTOMRIGHT:
		A.SubMatrix(ATL, 0, 0, k+mb, k+mb)
		A.SubMatrix(ATR, 0, k+mb, k+mb, A.Cols()-k-mb)

		A.SubMatrix(ABL, k+mb, 0, A.Rows()-k-mb, k+mb)
		A.SubMatrix(ABR, k+mb, k+mb)
	case pTOPLEFT:
		A.SubMatrix(ATL, 0, 0, k, k)
		A.SubMatrix(ATR, 0, k, k, A.Cols()-k)

		A.SubMatrix(ABL, k, 0, A.Rows()-k, A.Cols()-k)
		A.SubMatrix(ABR, k, k)
	}
}

/*
 * Partition p to 2 by 1 blocks.
 *
 *        pT
 *  p --> --
 *        pB
 *
 * Parameter nb is initial block size for pT (pTOP) or pB (pBOTTOM).
 */
func partitionPivot2x1(pT, pB, p *pPivots, nb int, pdir pDirection) {
	switch pdir {
	case pTOP:
		if nb == 0 {
			pT.pivots = nil
		} else {
			pT.pivots = p.pivots[:nb]
		}
		pB.pivots = p.pivots[nb:]
	case pBOTTOM:
		if nb > 0 {
			pT.pivots = p.pivots[:-nb]
			pT.pivots = p.pivots[len(p.pivots)-nb:]
		} else {
			pT.pivots = p.pivots
			pB.pivots = nil
		}
	}
}

/*
 * Repartition 2 by 1 block to 3 by 1 block.
 *
 *           pT      p0            pT       p0
 * pBOTTOM: --  --> --   ; pTOP:   --  -->  p1
 *           pB      p1            pB       --
 *                   p2                     p2
 *
 */
func repartPivot2x1to3x1(pT, p0, p1, p2, p *pPivots, nb int, pdir pDirection) {
	nT := len(pT.pivots)
	switch pdir {
	case pBOTTOM:
		if nT+nb > len(p.pivots) {
			nb = len(p.pivots) - nT
		}
		p0.pivots = pT.pivots
		p1.pivots = p.pivots[nT : nT+nb]
		p2.pivots = p.pivots[nT+nb:]
	case pTOP:
		if nb > nT {
			nb = nT
		}
		p0.pivots = p.pivots[:nT-nb]
		p1.pivots = p.pivots[nT-nb : nT]
		p2.pivots = p.pivots[nT:]
	}
}

/*
 * Continue with 2 by 1 block from 3 by 1 block.
 *
 *           pT      p0            pT       p0
 * pBOTTOM: --  <--  p1   ; pTOP:   -- <--  --
 *           pB      --            pB       p1
 *                   p2                     p2
 *
 */
func contPivot3x1to2x1(pT, pB, p0, p1, p *pPivots, pdir pDirection) {
	var n0, n1 int
	n0 = len(p0.pivots)
	n1 = len(p1.pivots)
	switch pdir {
	case pBOTTOM:
		pT.pivots = p.pivots[:n0+n1]
		pB.pivots = p.pivots[n0+n1:]
	case pTOP:
		pT.pivots = p.pivots[:n0]
		pB.pivots = p.pivots[n0:]
	}
}

// Local Variables:
// tab-width: 4
// indent-tabs-mode: nil
// End:
