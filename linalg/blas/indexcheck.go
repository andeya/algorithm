// Copyright (c) Harri Rautila, 2012,2013

// This file is part of github.com/hrautila/linalg/blas package.
// It is free software, distributed under the terms of GNU Lesser General Public
// License Version 3, or any later version. See the COPYING tile included in this archive.

package blas

import (
	"errors"
	"github.com/henrylee2cn/algorithm/linalg"
	"github.com/henrylee2cn/algorithm/matrix"
)

type funcNum int

const (
	fnrm2 = 1 + iota
	fasum
	fiamax
	fdot
	fswap
	fcopy
	fset
	faxpy
	faxpby
	fscal
	frotg
	frotmg
	frot
	frotm
	fgemv
	fgbmv
	fdtrm
	ftbmv
	ftpmv
	ftrsv
	ftbsv
	ftpsv
	ftrmv
	fsymv
	fsbmv
	fspmv
	fger
	fsyr
	fspr
	fsyr2
	fdspr2
	fgemm
	fsymm
	fsyrk
	fsyr2k
	ftrmm
	ftrsm
)

func abs(val int) int {
	if val < 0 {
		return -val
	}
	return val
}

func max(a, b int) int {
	if a < b {
		return b
	}
	return a
}

var panicOnError bool = false

func PanicOnError(flag bool) {
	panicOnError = flag
}

func onError(msg string) error {
	if panicOnError {
		panic(msg)
	}
	return errors.New(msg)
}

func check_level1_func(ind *linalg.IndexOpts, fn funcNum, X, Y matrix.Matrix) error {

	nX, nY := 0, 0
	// this is adapted from cvxopt:blas.c python blas interface
	switch fn {
	case fnrm2, fasum, fiamax, fscal, fset:
		if ind.IncX <= 0 {
			return onError("incX illegal, <=0")
		}
		if ind.OffsetX < 0 {
			return onError("offsetX illegal, <0")
		}
		sizeX := X.NumElements()
		if sizeX >= ind.OffsetX+1 {
			// calculate default size for N based on X size
			nX = 1 + (sizeX-ind.OffsetX-1)/ind.IncX
		}
		if sizeX < ind.OffsetX+1+(ind.Nx-1)*abs(ind.IncX) {
			return onError("X size error")
		}
		if ind.Nx < 0 {
			ind.Nx = nX
		}

	case fdot, fswap, fcopy, faxpy, faxpby:
		// vector X
		if ind.IncX <= 0 {
			return onError("incX illegal, <=0")
		}
		if ind.OffsetX < 0 {
			return onError("offsetX illegal, <0")
		}
		sizeX := X.NumElements()
		if sizeX >= ind.OffsetX+1 {
			// calculate default size for N based on X size
			nX = 1 + (sizeX-ind.OffsetX-1)/ind.IncX
		}
		if sizeX < ind.OffsetX+1+(ind.Nx-1)*abs(ind.IncX) {
			return onError("X size error")
		}
		if ind.Nx < 0 {
			ind.Nx = nX
		}
		// vector Y
		if ind.IncY <= 0 {
			return onError("incY illegal, <=0")
		}
		if ind.OffsetY < 0 {
			return onError("offsetY illegal, <0")
		}
		sizeY := Y.NumElements()
		if sizeY >= ind.OffsetY+1 {
			// calculate default size for N based on Y size
			nY = 1 + (sizeY-ind.OffsetY-1)/ind.IncY
		}
		if ind.Ny < 0 {
			ind.Ny = nY
		}
		if sizeY < ind.OffsetY+1+(ind.Ny-1)*abs(ind.IncY) {
			//fmt.Printf("sizeY=%d, inds: %#v\n", sizeY, ind)
			return onError("Y size error")
		}

	case frotg, frotmg, frot, frotm:
	}
	return nil
}

func check_level2_func(ind *linalg.IndexOpts, fn funcNum, X, Y, A matrix.Matrix, pars *linalg.Parameters) error {
	if ind.IncX <= 0 {
		return onError("incX")
	}
	if ind.IncY <= 0 {
		return onError("incY")
	}

	sizeA := A.NumElements()
	arows := ind.LDa
	switch fn {
	case fgemv: // general matrix
		if ind.M < 0 {
			ind.M = A.Rows()
		}
		if ind.N < 0 {
			ind.N = A.Cols()
		}
		if ind.LDa == 0 {
			ind.LDa = max(1, A.LeadingIndex())
			arows = max(1, A.Rows())
		}
		if ind.OffsetA < 0 {
			return onError("offsetA")
		}
		if ind.N > 0 && ind.M > 0 &&
			sizeA < ind.OffsetA+(ind.N-1)*arows+ind.M {
			return onError("sizeA")
		}
		if ind.OffsetX < 0 {
			return onError("offsetX")
		}
		if ind.OffsetY < 0 {
			return onError("offsetY")
		}
		sizeX := X.NumElements()
		sizeY := Y.NumElements()
		if pars.Trans == linalg.PNoTrans {
			if ind.N > 0 && sizeX < ind.OffsetX+(ind.N-1)*abs(ind.IncX)+1 {
				return onError("sizeX")
			}
			if ind.M > 0 && sizeY < ind.OffsetY+(ind.M-1)*abs(ind.IncY)+1 {
				return onError("sizeY")
			}
		} else {
			if ind.M > 0 && sizeX < ind.OffsetX+(ind.M-1)*abs(ind.IncX)+1 {
				return onError("sizeX")
			}
			if ind.N > 0 && sizeY < ind.OffsetY+(ind.N-1)*abs(ind.IncY)+1 {
				return onError("sizeY")
			}
		}
	case fger:
		if ind.M < 0 {
			ind.M = A.Rows()
		}
		if ind.N < 0 {
			ind.N = A.Cols()
		}
		if ind.M == 0 || ind.N == 0 {
			return nil
		}
		if ind.M > 0 && ind.N > 0 {
			if ind.LDa == 0 {
				ind.LDa = max(1, A.LeadingIndex())
				arows = max(1, A.Rows())
			}
			if ind.LDa < max(1, ind.M) {
				return onError("ldA")
			}
			if ind.OffsetA < 0 {
				return onError("offsetA")
			}
			if sizeA < ind.OffsetA+(ind.N-1)*arows+ind.M {
				return onError("sizeA")
			}
			if ind.OffsetX < 0 {
				return onError("offsetX")
			}
			if ind.OffsetY < 0 {
				return onError("offsetY")
			}
			sizeX := X.NumElements()
			if sizeX < ind.OffsetX+(ind.M-1)*abs(ind.IncX)+1 {
				return onError("sizeX")
			}
			sizeY := Y.NumElements()
			if sizeY < ind.OffsetY+(ind.N-1)*abs(ind.IncY)+1 {
				return onError("sizeY")
			}
		}
	case fgbmv: // general banded
		if ind.M < 0 {
			ind.M = A.Rows()
		}
		if ind.N < 0 {
			ind.N = A.Cols()
		}
		if ind.Kl < 0 {
			return onError("kl")
		}
		if ind.Ku < 0 {
			ind.Ku = A.Rows() - 1 - ind.Kl
		}
		if ind.Ku < 0 {
			return onError("ku")
		}
		if ind.LDa == 0 {
			ind.LDa = max(1, A.LeadingIndex())
			arows = max(1, A.Rows())
		}
		if ind.LDa < ind.Kl+ind.Ku+1 {
			return onError("ldA")
		}
		if ind.OffsetA < 0 {
			return onError("offsetA")
		}
		sizeA := A.NumElements()
		if ind.N > 0 && ind.M > 0 &&
			sizeA < ind.OffsetA+(ind.N-1)*arows+ind.Kl+ind.Ku+1 {
			return onError("sizeA")
		}
		if ind.OffsetX < 0 {
			return onError("offsetX")
		}
		if ind.OffsetY < 0 {
			return onError("offsetY")
		}
		sizeX := X.NumElements()
		sizeY := Y.NumElements()
		if pars.Trans == linalg.PNoTrans {
			if ind.N > 0 && sizeX < ind.OffsetX+(ind.N-1)*abs(ind.IncX)+1 {
				return onError("sizeX")
			}
			if ind.N > 0 && sizeY < ind.OffsetY+(ind.M-1)*abs(ind.IncY)+1 {
				return onError("sizeY")
			}
		} else {
			if ind.N > 0 && sizeX < ind.OffsetX+(ind.M-1)*abs(ind.IncX)+1 {
				return onError("sizeX")
			}
			if ind.N > 0 && sizeY < ind.OffsetY+(ind.N-1)*abs(ind.IncY)+1 {
				return onError("sizeY")
			}
		}
	case ftrmv, ftrsv:
		// ftrmv = triangular
		// ftrsv = triangular solve
		if ind.N < 0 {
			if A.Rows() != A.Cols() {
				return onError("A not square")
			}
			ind.N = A.Rows()
		}
		if ind.N > 0 {
			if ind.LDa == 0 {
				ind.LDa = max(1, A.LeadingIndex())
				arows = max(1, A.Rows())
			}
			if ind.LDa < max(1, ind.N) {
				return onError("ldA")
			}
			if ind.OffsetA < 0 {
				return onError("offsetA")
			}
			sizeA := A.NumElements()
			if sizeA < ind.OffsetA+(ind.N-1)*arows+ind.N {
				return onError("sizeA")
			}
			sizeX := X.NumElements()
			if sizeX < ind.OffsetX+(ind.N-1)*abs(ind.IncX)+1 {
				return onError("sizeX")
			}
		}
	case ftbmv, ftbsv, fsbmv:
		// ftbmv = triangular banded
		// ftbsv = triangular banded solve
		// fsbmv = symmetric banded product
		arows := ind.LDa
		if ind.N < 0 {
			ind.N = A.Rows()
		}
		if ind.N > 0 {
			if ind.K < 0 {
				ind.K = max(0, A.Rows()-1)
			}
			if ind.LDa == 0 {
				ind.LDa = max(1, A.LeadingIndex())
				arows = max(1, A.Rows())
			}
			if ind.LDa < ind.K+1 {
				return onError("ldA")
			}
			if ind.OffsetA < 0 {
				return onError("offsetA")
			}
			sizeA := A.NumElements()
			if sizeA < ind.OffsetA+(ind.N-1)*arows+ind.K+1 {
				return onError("sizeA")
			}
			sizeX := X.NumElements()
			if sizeX < ind.OffsetX+(ind.N-1)*abs(ind.IncX)+1 {
				return onError("sizeX")
			}
			if Y != nil {
				sizeY := Y.NumElements()
				if sizeY < ind.OffsetY+(ind.N-1)*abs(ind.IncY)+1 {
					return onError("sizeY")
				}
			}
		}
	case fsymv, fsyr, fsyr2:
		// fsymv = symmetric product
		// fsyr = symmetric rank update
		// fsyr2 = symmetric rank-2 update
		if ind.N < 0 {
			if A.Rows() != A.Cols() {
				return onError("A not square")
			}
			ind.N = A.Rows()
		}
		arows := ind.LDa
		if ind.N > 0 {
			if ind.LDa == 0 {
				ind.LDa = max(1, A.LeadingIndex())
				arows = max(1, A.Rows())
			}
			if ind.LDa < max(1, ind.N) {
				return onError("ldA")
			}
			if ind.OffsetA < 0 {
				return onError("offsetA")
			}
			sizeA := A.NumElements()
			if sizeA < ind.OffsetA+(ind.N-1)*arows+ind.N {
				return onError("sizeA")
			}
			if ind.OffsetX < 0 {
				return onError("offsetX")
			}
			sizeX := X.NumElements()
			if sizeX < ind.OffsetX+(ind.N-1)*abs(ind.IncX)+1 {
				return onError("sizeX")
			}
			if Y != nil {
				if ind.OffsetY < 0 {
					return onError("offsetY")
				}
				sizeY := Y.NumElements()
				if sizeY < ind.OffsetY+(ind.N-1)*abs(ind.IncY)+1 {
					return onError("sizeY")
				}
			}
		}
	case fspr, fdspr2, ftpsv, fspmv, ftpmv:
		// ftpsv = triangular packed solve
		// fspmv = symmetric packed product
		// ftpmv = triangular packed
	}
	return nil
}

func check_level3_func(ind *linalg.IndexOpts, fn funcNum, A, B, C matrix.Matrix,
	pars *linalg.Parameters) (err error) {

	// defaults for these
	arows := ind.LDa
	brows := ind.LDb
	crows := ind.LDc

	switch fn {
	case fgemm:
		if ind.M < 0 {
			if pars.TransA == linalg.PNoTrans {
				ind.M = A.Rows()
			} else {
				ind.M = A.Cols()
			}
		}
		if ind.N < 0 {
			if pars.TransB == linalg.PNoTrans {
				ind.N = B.Cols()
			} else {
				ind.N = B.Rows()
			}
		}
		if ind.M == 0 || ind.N == 0 {
			return nil
		}
		if ind.K < 0 {
			if pars.TransA == linalg.PNoTrans {
				ind.K = A.Cols()
			} else {
				ind.K = A.Rows()
			}
			if pars.TransB == linalg.PNoTrans && ind.K != B.Rows() ||
				pars.TransB != linalg.PNoTrans && ind.K != B.Cols() {
				return onError("dimensions of A and B do not match")
			}
		}
		if ind.OffsetA < 0 {
			return onError("offsetA illegal, <0")
		}
		if ind.LDa == 0 {
			ind.LDa = max(1, A.LeadingIndex())
			arows = max(1, A.Rows())
		}
		if ind.K > 0 {
			if (pars.TransA == linalg.PNoTrans && ind.LDa < max(1, ind.M)) ||
				(pars.TransA != linalg.PNoTrans && ind.LDa < max(1, ind.K)) {
				return onError("inconsistent ldA")
			}
			sizeA := A.NumElements()
			if (pars.TransA == linalg.PNoTrans &&
				sizeA < ind.OffsetA+(ind.K-1)*arows+ind.M) ||
				(pars.TransA != linalg.PNoTrans &&
					sizeA < ind.OffsetA+(ind.M-1)*arows+ind.K) {
				return onError("sizeA")
			}
		}
		// B matrix
		if ind.OffsetB < 0 {
			return onError("offsetB illegal, <0")
		}
		if ind.LDb == 0 {
			ind.LDb = max(1, B.LeadingIndex())
			brows = max(1, B.Rows())
		}
		if ind.K > 0 {
			if (pars.TransB == linalg.PNoTrans && ind.LDb < max(1, ind.K)) ||
				(pars.TransB != linalg.PNoTrans && ind.LDb < max(1, ind.N)) {
				return onError("inconsistent ldB")
			}
			sizeB := B.NumElements()
			if (pars.TransB == linalg.PNoTrans &&
				sizeB < ind.OffsetB+(ind.N-1)*brows+ind.K) ||
				(pars.TransB != linalg.PNoTrans &&
					sizeB < ind.OffsetB+(ind.K-1)*brows+ind.N) {
				return onError("sizeB")
			}
		}
		// C matrix
		if ind.OffsetC < 0 {
			return onError("offsetC illegal, <0")
		}
		if ind.LDc == 0 {
			ind.LDc = max(1, C.LeadingIndex())
			crows = max(1, C.Rows())
		}
		if ind.LDc < max(1, ind.M) {
			return onError("inconsistent ldC")
		}
		sizeC := C.NumElements()
		if sizeC < ind.OffsetC+(ind.N-1)*crows+ind.M {
			return onError("sizeC")
		}

	case fsymm, ftrmm, ftrsm:
		if ind.M < 0 {
			ind.M = B.Rows()
			if pars.Side == linalg.PLeft && (ind.M != A.Rows() || ind.M != A.Cols()) {
				return onError("dimensions of A and B do not match")
			}
		}
		if ind.N < 0 {
			ind.N = B.Cols()
			if pars.Side == linalg.PRight && (ind.N != A.Rows() || ind.N != A.Cols()) {
				return onError("dimensions of A and B do not match")
			}
		}
		if ind.M == 0 || ind.N == 0 {
			return
		}
		// check A
		if ind.OffsetB < 0 {
			return onError("offsetB illegal, <0")
		}
		if ind.LDa == 0 {
			ind.LDa = max(1, A.LeadingIndex())
			arows = max(1, A.Rows())
		}
		if pars.Side == linalg.PLeft && ind.LDa < max(1, ind.M) || ind.LDa < max(1, ind.N) {
			return onError("ldA")
		}
		sizeA := A.NumElements()
		if (pars.Side == linalg.PLeft && sizeA < ind.OffsetA+(ind.M-1)*arows+ind.M) ||
			(pars.Side == linalg.PRight && sizeA < ind.OffsetA+(ind.N-1)*arows+ind.N) {
			return onError("sizeA")
		}

		if B != nil {
			if ind.OffsetB < 0 {
				return onError("offsetB illegal, <0")
			}
			if ind.LDb == 0 {
				ind.LDb = max(1, B.LeadingIndex())
				brows = max(1, B.Rows())
			}
			if ind.LDb < max(1, ind.M) {
				return onError("ldB")
			}
			sizeB := B.NumElements()
			if sizeB < ind.OffsetB+(ind.N-1)*brows+ind.M {
				return onError("sizeB")
			}
		}

		if C != nil {
			if ind.OffsetC < 0 {
				return onError("offsetC illegal, <0")
			}
			if ind.LDc == 0 {
				ind.LDc = max(1, C.LeadingIndex())
				crows = max(1, C.Rows())
			}
			if ind.LDc < max(1, ind.M) {
				return onError("ldC")
			}
			sizeC := C.NumElements()
			if sizeC < ind.OffsetC+(ind.N-1)*crows+ind.M {
				return onError("sizeC")
			}
		}
	case fsyrk:
		if ind.N < 0 {
			if pars.Trans == linalg.PNoTrans {
				ind.N = A.Rows()
			} else {
				ind.N = A.Cols()
			}
			//ind.N = C.Rows()
		}
		if ind.K < 0 {
			if pars.Trans == linalg.PNoTrans {
				ind.K = A.Cols()
			} else {
				ind.K = A.Rows()
			}
		}
		if ind.N == 0 {
			return
		}
		if ind.LDa == 0 {
			ind.LDa = max(1, A.LeadingIndex())
			arows = max(1, A.Rows())
		}
		if ind.OffsetA < 0 {
			return onError("offsetA")
		}
		if ind.K > 0 {
			if (pars.Trans == linalg.PNoTrans && ind.LDa < max(1, ind.N)) ||
				(pars.Trans != linalg.PNoTrans && ind.LDa < max(1, ind.K)) {
				return onError("inconsistent ldA")
			}
			sizeA := A.NumElements()
			if (pars.Trans == linalg.PNoTrans &&
				sizeA < ind.OffsetA+(ind.K-1)*arows+ind.N) ||
				(pars.TransA != linalg.PNoTrans &&
					sizeA < ind.OffsetA+(ind.N-1)*arows+ind.K) {
				return onError("sizeA")
			}
		}

		if ind.OffsetC < 0 {
			return onError("offsetC illegal, <0")
		}
		if ind.LDc == 0 {
			ind.LDc = max(1, C.LeadingIndex())
			crows = max(1, C.Rows())
		}
		if ind.LDc < max(1, ind.N) {
			return onError("ldC")
		}
		sizeC := C.NumElements()
		if sizeC < ind.OffsetC+(ind.N-1)*crows+ind.N {
			return onError("sizeC")
		}
	case fsyr2k:
		if ind.N < 0 {
			if pars.Trans == linalg.PNoTrans {
				ind.N = A.Rows()
				if ind.N != B.Rows() {
					return onError("dimensions of A and B do not match")
				}
			} else {
				ind.N = A.Cols()
				if ind.N != B.Cols() {
					return onError("dimensions of A and B do not match")
				}
			}
		}
		if ind.N == 0 {
			return
		}
		if ind.K < 0 {
			if pars.Trans == linalg.PNoTrans {
				ind.K = A.Cols()
				if ind.K != B.Cols() {
					return onError("dimensions of A and B do not match")
				}
			} else {
				ind.K = A.Rows()
				if ind.K != B.Rows() {
					return onError("dimensions of A and B do not match")
				}
			}
		}
		if ind.LDa == 0 {
			ind.LDa = max(1, A.LeadingIndex())
			arows = max(1, A.Rows())
		}
		if ind.K > 0 {
			if (pars.Trans == linalg.PNoTrans && ind.LDa < max(1, ind.N)) ||
				(pars.Trans != linalg.PNoTrans && ind.LDa < max(1, ind.K)) {
				return onError("inconsistent ldA")
			}
			sizeA := A.NumElements()
			if (pars.Trans == linalg.PNoTrans &&
				sizeA < ind.OffsetA+(ind.K-1)*arows+ind.N) ||
				(pars.TransA != linalg.PNoTrans &&
					sizeA < ind.OffsetA+(ind.N-1)*arows+ind.K) {
				return onError("sizeA")
			}
		}
		if ind.OffsetB < 0 {
			return onError("offsetB illegal, <0")
		}
		if ind.LDb == 0 {
			ind.LDb = max(1, B.LeadingIndex())
			brows = max(1, B.Rows())
		}
		if ind.K > 0 {
			if (pars.Trans == linalg.PNoTrans && ind.LDb < max(1, ind.N)) ||
				(pars.Trans != linalg.PNoTrans && ind.LDb < max(1, ind.K)) {
				return onError("ldB")
			}
			sizeB := B.NumElements()
			if (pars.Trans == linalg.PNoTrans &&
				sizeB < ind.OffsetB+(ind.K-1)*brows+ind.N) ||
				(pars.Trans != linalg.PNoTrans &&
					sizeB < ind.OffsetB+(ind.N-1)*brows+ind.K) {
				return onError("sizeB")
			}
		}
		if ind.OffsetC < 0 {
			return onError("offsetC illegal, <0")
		}
		if ind.LDc == 0 {
			ind.LDc = max(1, C.LeadingIndex())
			crows = max(1, C.Rows())
		}
		if ind.LDc < max(1, ind.N) {
			return onError("ldC")
		}
		sizeC := C.NumElements()
		if sizeC < ind.OffsetC+(ind.N-1)*crows+ind.N {
			return onError("sizeC")
		}
	}
	err = nil
	return
}

// Local Variables:
// tab-width: 4
// End:
