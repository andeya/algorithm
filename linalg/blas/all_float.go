// Copyright (c) Harri Rautila, 2012,2013

// This file is part of github.com/hrautila/linalg package.
// It is free software, distributed under the terms of GNU Lesser General Public
// License Version 3, or any later version. See the COPYING tile included in this archive.

package blas

import (
	"github.com/henrylee2cn/algorithm/linalg"
	"github.com/henrylee2cn/algorithm/matrix"
	"math"
)

// See function Nrm2.
func Nrm2Float(X *matrix.FloatMatrix, opts ...linalg.Option) (v float64) {
	v = math.NaN()
	ind := linalg.GetIndexOpts(opts...)
	err := check_level1_func(ind, fnrm2, X, nil)
	if err != nil {
		return
	}
	if ind.Nx == 0 {
		v = 0.0
		return
	}
	Xa := X.FloatArray()
	v = dnrm2(ind.Nx, Xa[ind.OffsetX:], ind.IncX)
	return
}

// See function Asum.
func AsumFloat(X *matrix.FloatMatrix, opts ...linalg.Option) (v float64) {
	v = math.NaN()
	ind := linalg.GetIndexOpts(opts...)
	err := check_level1_func(ind, fasum, X, nil)
	if err != nil {
		return
	}
	if ind.Nx == 0 {
		v = 0.0
		return
	}
	Xa := X.FloatArray()
	v = dasum(ind.Nx, Xa[ind.OffsetX:], ind.IncX)
	return
}

// See functin Dot.
func DotFloat(X, Y *matrix.FloatMatrix, opts ...linalg.Option) (v float64) {
	v = math.NaN()
	ind := linalg.GetIndexOpts(opts...)
	err := check_level1_func(ind, fdot, X, Y)
	if err != nil {
		return
	}
	if ind.Nx == 0 {
		v = 0.0
		return
	}
	Xa := X.FloatArray()
	Ya := Y.FloatArray()
	v = ddot(ind.Nx, Xa[ind.OffsetX:], ind.IncX, Ya[ind.OffsetY:], ind.IncY)
	return
}

// See function Swap.
func SwapFloat(X, Y *matrix.FloatMatrix, opts ...linalg.Option) (err error) {
	ind := linalg.GetIndexOpts(opts...)
	err = check_level1_func(ind, fswap, X, Y)
	if err != nil {
		return
	}
	if ind.Nx == 0 {
		return
	}
	Xa := X.FloatArray()
	Ya := Y.FloatArray()
	dswap(ind.Nx, Xa[ind.OffsetX:], ind.IncX, Ya[ind.OffsetY:], ind.IncY)
	return
}

// See function Copy.
func CopyFloat(X, Y *matrix.FloatMatrix, opts ...linalg.Option) (err error) {
	ind := linalg.GetIndexOpts(opts...)
	err = check_level1_func(ind, fcopy, X, Y)
	if err != nil {
		return
	}
	if ind.Nx == 0 {
		return
	}
	Xa := X.FloatArray()
	Ya := Y.FloatArray()
	dcopy(ind.Nx, Xa[ind.OffsetX:], ind.IncX, Ya[ind.OffsetY:], ind.IncY)
	return
}

// Calculate for column vector X := alpha * X. Contents of X is set to new value.
// Valid options: offset, inc, N.

// See function Scal.
func ScalFloat(X *matrix.FloatMatrix, alpha float64, opts ...linalg.Option) (err error) {
	ind := linalg.GetIndexOpts(opts...)
	err = check_level1_func(ind, fscal, X, nil)
	if err != nil {
		return
	}
	if ind.Nx == 0 {
		return
	}
	Xa := X.FloatArray()
	dscal(ind.Nx, alpha, Xa[ind.OffsetX:], ind.IncX)
	return
}

// See function Axpy.
func AxpyFloat(X, Y *matrix.FloatMatrix, alpha float64, opts ...linalg.Option) (err error) {
	ind := linalg.GetIndexOpts(opts...)
	err = check_level1_func(ind, faxpy, X, Y)
	if err != nil {
		return
	}
	if ind.Nx == 0 {
		return
	}
	Xa := X.FloatArray()
	Ya := Y.FloatArray()
	daxpy(ind.Nx, alpha, Xa[ind.OffsetX:], ind.IncX, Ya[ind.OffsetY:], ind.IncY)
	return
}

// ---------------------------------------------------------------------------------
// BLAS LEVEL 2
// ---------------------------------------------------------------------------------

// See function Gemv.
func GemvFloat(A, X, Y *matrix.FloatMatrix, alpha, beta float64, opts ...linalg.Option) (err error) {

	var params *linalg.Parameters
	params, err = linalg.GetParameters(opts...)
	if err != nil {
		return
	}
	ind := linalg.GetIndexOpts(opts...)
	err = check_level2_func(ind, fgemv, X, Y, A, params)
	if err != nil {
		return
	}
	if ind.M == 0 && params.Trans == linalg.PNoTrans {
		return
	}
	if ind.N == 0 && (params.Trans == linalg.PTrans || params.Trans == linalg.PConjTrans) {
		return
	}
	Xa := X.FloatArray()
	Ya := Y.FloatArray()
	Aa := A.FloatArray()
	if params.Trans == linalg.PNoTrans && ind.N == 0 {
		dscal(ind.M, beta, Ya[ind.OffsetY:], ind.IncY)
	} else if params.Trans == linalg.PTrans && ind.M == 0 {
		dscal(ind.N, beta, Ya[ind.OffsetY:], ind.IncY)
	} else {
		trans := linalg.ParamString(params.Trans)
		dgemv(trans, ind.M, ind.N, alpha, Aa[ind.OffsetA:],
			ind.LDa, Xa[ind.OffsetX:], ind.IncX, beta, Ya[ind.OffsetY:], ind.IncY)
	}
	return
}

// See function Gbmv.
func GbmvFloat(A, X, Y *matrix.FloatMatrix, alpha, beta float64, opts ...linalg.Option) (err error) {

	var params *linalg.Parameters
	params, err = linalg.GetParameters(opts...)
	if err != nil {
		return
	}
	ind := linalg.GetIndexOpts(opts...)
	err = check_level2_func(ind, fgbmv, X, Y, A, params)
	if err != nil {
		return
	}
	if ind.M == 0 && ind.N == 0 {
		return
	}

	Xa := X.FloatArray()
	Ya := Y.FloatArray()
	Aa := A.FloatArray()
	if params.Trans == linalg.PNoTrans && ind.N == 0 {
		dscal(ind.M, beta, Ya[ind.OffsetY:], ind.IncY)
	} else if params.Trans == linalg.PTrans && ind.M == 0 {
		dscal(ind.N, beta, Ya[ind.OffsetY:], ind.IncY)
	} else {
		trans := linalg.ParamString(params.Trans)
		dgbmv(trans, ind.M, ind.N, ind.Kl, ind.Ku,
			alpha, Aa[ind.OffsetA:], ind.LDa, Xa[ind.OffsetX:], ind.IncX,
			beta, Ya[ind.OffsetY:], ind.IncY)
	}
	return
}

// See function Symv.
func SymvFloat(A, X, Y *matrix.FloatMatrix, alpha, beta float64, opts ...linalg.Option) (err error) {

	var params *linalg.Parameters
	params, err = linalg.GetParameters(opts...)
	if err != nil {
		return
	}
	ind := linalg.GetIndexOpts(opts...)
	err = check_level2_func(ind, fsymv, X, Y, A, params)
	if err != nil {
		return
	}
	if ind.N == 0 {
		return
	}
	Xa := X.FloatArray()
	Ya := Y.FloatArray()
	Aa := A.FloatArray()
	uplo := linalg.ParamString(params.Uplo)
	dsymv(uplo, ind.N, alpha, Aa[ind.OffsetA:], ind.LDa, Xa[ind.OffsetX:], ind.IncX,
		beta, Ya[ind.OffsetY:], ind.IncY)
	return
}

// See function Sbmv.
func SbmvFloat(A, X, Y *matrix.FloatMatrix, alpha, beta float64, opts ...linalg.Option) (err error) {

	var params *linalg.Parameters
	params, err = linalg.GetParameters(opts...)
	if err != nil {
		return
	}
	ind := linalg.GetIndexOpts(opts...)
	err = check_level2_func(ind, fsbmv, X, Y, A, params)
	if err != nil {
		return
	}
	if ind.N == 0 {
		return
	}
	Xa := X.FloatArray()
	Ya := Y.FloatArray()
	Aa := A.FloatArray()
	uplo := linalg.ParamString(params.Uplo)
	dsbmv(uplo, ind.N, ind.K, alpha, Aa[ind.OffsetA:], ind.LDa, Xa[ind.OffsetX:],
		ind.IncX, beta, Ya[ind.OffsetY:], ind.IncY)
	return
}

// See function Trmv.
func TrmvFloat(A, X *matrix.FloatMatrix, opts ...linalg.Option) (err error) {

	var params *linalg.Parameters
	params, err = linalg.GetParameters(opts...)
	if err != nil {
		return
	}
	ind := linalg.GetIndexOpts(opts...)
	err = check_level2_func(ind, ftrmv, X, nil, A, params)
	if err != nil {
		return
	}
	if ind.N == 0 {
		return
	}
	Xa := X.FloatArray()
	Aa := A.FloatArray()
	uplo := linalg.ParamString(params.Uplo)
	trans := linalg.ParamString(params.Trans)
	diag := linalg.ParamString(params.Diag)
	dtrmv(uplo, trans, diag, ind.N,
		Aa[ind.OffsetA:], ind.LDa, Xa[ind.OffsetX:], ind.IncX)
	return
}

// See function Tbmv.
func TbmvFloat(A, X *matrix.FloatMatrix, opts ...linalg.Option) (err error) {

	var params *linalg.Parameters
	params, err = linalg.GetParameters(opts...)
	if err != nil {
		return
	}
	ind := linalg.GetIndexOpts(opts...)
	err = check_level2_func(ind, ftbmv, X, nil, A, params)
	if err != nil {
		return
	}
	if ind.N == 0 {
		return
	}
	Xa := X.FloatArray()
	Aa := A.FloatArray()
	uplo := linalg.ParamString(params.Uplo)
	trans := linalg.ParamString(params.Trans)
	diag := linalg.ParamString(params.Diag)
	dtbmv(uplo, trans, diag, ind.N, ind.K,
		Aa[ind.OffsetA:], ind.LDa, Xa[ind.OffsetX:], ind.IncX)
	return
}

// See function Trsv.
func TrsvFloat(A, X *matrix.FloatMatrix, opts ...linalg.Option) (err error) {

	var params *linalg.Parameters
	params, err = linalg.GetParameters(opts...)
	if err != nil {
		return
	}
	ind := linalg.GetIndexOpts(opts...)
	err = check_level2_func(ind, ftrsv, X, nil, A, params)
	if err != nil {
		return
	}
	if ind.N == 0 {
		return
	}
	Xa := X.FloatArray()
	Aa := A.FloatArray()
	uplo := linalg.ParamString(params.Uplo)
	trans := linalg.ParamString(params.Trans)
	diag := linalg.ParamString(params.Diag)
	dtrsv(uplo, trans, diag, ind.N,
		Aa[ind.OffsetA:], ind.LDa, Xa[ind.OffsetX:], ind.IncX)
	return
}

// See function Tbsv.
func TbsvFloat(A, X *matrix.FloatMatrix, opts ...linalg.Option) (err error) {

	var params *linalg.Parameters
	params, err = linalg.GetParameters(opts...)
	if err != nil {
		return
	}
	ind := linalg.GetIndexOpts(opts...)
	err = check_level2_func(ind, ftbsv, X, nil, A, params)
	if err != nil {
		return
	}
	if ind.N == 0 {
		return
	}
	Xa := X.FloatArray()
	Aa := A.FloatArray()
	uplo := linalg.ParamString(params.Uplo)
	trans := linalg.ParamString(params.Trans)
	diag := linalg.ParamString(params.Diag)
	dtbsv(uplo, trans, diag, ind.N, ind.K,
		Aa[ind.OffsetA:], ind.LDa, Xa[ind.OffsetX:], ind.IncX)
	return
}

// See function Ger.
func GerFloat(X, Y, A *matrix.FloatMatrix, alpha float64, opts ...linalg.Option) (err error) {

	var params *linalg.Parameters
	params, err = linalg.GetParameters(opts...)
	if err != nil {
		return
	}
	ind := linalg.GetIndexOpts(opts...)
	err = check_level2_func(ind, fger, X, Y, A, params)
	if err != nil {
		return
	}
	if ind.N == 0 || ind.M == 0 {
		return
	}
	Xa := X.FloatArray()
	Ya := Y.FloatArray()
	Aa := A.FloatArray()
	dger(ind.M, ind.N, alpha, Xa[ind.OffsetX:], ind.IncX,
		Ya[ind.OffsetY:], ind.IncY, Aa[ind.OffsetA:], ind.LDa)

	return
}

// See function Syr.
func SyrFloat(X, A *matrix.FloatMatrix, alpha float64, opts ...linalg.Option) (err error) {

	var params *linalg.Parameters
	params, err = linalg.GetParameters(opts...)
	if err != nil {
		return
	}
	ind := linalg.GetIndexOpts(opts...)
	err = check_level2_func(ind, fsyr, X, nil, A, params)
	if err != nil {
		return
	}
	if ind.N == 0 {
		return
	}
	Xa := X.FloatArray()
	Aa := A.FloatArray()
	uplo := linalg.ParamString(params.Uplo)
	dsyr(uplo, ind.N, alpha, Xa[ind.OffsetX:], ind.IncX, Aa[ind.OffsetA:], ind.LDa)
	return
}

// See function Syr2.
func Syr2Float(X, Y, A *matrix.FloatMatrix, alpha float64, opts ...linalg.Option) (err error) {

	var params *linalg.Parameters
	params, err = linalg.GetParameters(opts...)
	if err != nil {
		return
	}
	ind := linalg.GetIndexOpts(opts...)
	err = check_level2_func(ind, fsyr2, X, Y, A, params)
	if err != nil {
		return
	}
	if ind.N == 0 {
		return
	}
	Xa := X.FloatArray()
	Ya := Y.FloatArray()
	Aa := A.FloatArray()
	uplo := linalg.ParamString(params.Uplo)
	dsyr2(uplo, ind.N, alpha, Xa[ind.OffsetX:], ind.IncX, Ya[ind.OffsetY:], ind.IncY,
		Aa[ind.OffsetA:], ind.LDa)
	return
}

// ---------------------------------------------------------------------------------
// BLAS LEVEL 3
// ---------------------------------------------------------------------------------

// See function Gemm.
func GemmFloat(A, B, C *matrix.FloatMatrix, alpha, beta float64, opts ...linalg.Option) (err error) {

	params, e := linalg.GetParameters(opts...)
	if e != nil {
		err = e
		return
	}
	ind := linalg.GetIndexOpts(opts...)
	err = check_level3_func(ind, fgemm, A, B, C, params)
	if err != nil {
		return
	}
	if ind.M == 0 || ind.N == 0 {
		return
	}
	Aa := A.FloatArray()
	Ba := B.FloatArray()
	Ca := C.FloatArray()
	transB := linalg.ParamString(params.TransB)
	transA := linalg.ParamString(params.TransA)
	//diag := linalg.ParamString(params.Diag)
	dgemm(transA, transB, ind.M, ind.N, ind.K, alpha,
		Aa[ind.OffsetA:], ind.LDa, Ba[ind.OffsetB:], ind.LDb, beta,
		Ca[ind.OffsetC:], ind.LDc)
	return
}

// See function Symm.
func SymmFloat(A, B, C *matrix.FloatMatrix, alpha, beta float64, opts ...linalg.Option) (err error) {

	params, e := linalg.GetParameters(opts...)
	if e != nil {
		err = e
		return
	}
	ind := linalg.GetIndexOpts(opts...)
	err = check_level3_func(ind, fsymm, A, B, C, params)
	if err != nil {
		return
	}
	if ind.M == 0 || ind.N == 0 {
		return
	}
	Aa := A.FloatArray()
	Ba := B.FloatArray()
	Ca := C.FloatArray()
	uplo := linalg.ParamString(params.Uplo)
	side := linalg.ParamString(params.Side)
	dsymm(side, uplo, ind.M, ind.N, alpha, Aa[ind.OffsetA:], ind.LDa,
		Ba[ind.OffsetB:], ind.LDb, beta, Ca[ind.OffsetC:], ind.LDc)

	return
}

// See function Syrk.
func SyrkFloat(A, C *matrix.FloatMatrix, alpha, beta float64, opts ...linalg.Option) (err error) {

	params, e := linalg.GetParameters(opts...)
	if e != nil {
		err = e
		return
	}
	ind := linalg.GetIndexOpts(opts...)
	err = check_level3_func(ind, fsyrk, A, nil, C, params)
	if e != nil || err != nil {
		return
	}
	if ind.N == 0 {
		return
	}
	Aa := A.FloatArray()
	Ca := C.FloatArray()
	uplo := linalg.ParamString(params.Uplo)
	trans := linalg.ParamString(params.Trans)
	//diag := linalg.ParamString(params.Diag)
	dsyrk(uplo, trans, ind.N, ind.K, alpha, Aa[ind.OffsetA:], ind.LDa, beta,
		Ca[ind.OffsetC:], ind.LDc)

	return
}

// See function Syrk2.
func Syr2kFloat(A, B, C *matrix.FloatMatrix, alpha, beta float64, opts ...linalg.Option) (err error) {

	params, e := linalg.GetParameters(opts...)
	if e != nil {
		err = e
		return
	}
	ind := linalg.GetIndexOpts(opts...)
	err = check_level3_func(ind, fsyr2k, A, B, C, params)
	if err != nil {
		return
	}
	if ind.N == 0 {
		return
	}
	Aa := A.FloatArray()
	Ba := B.FloatArray()
	Ca := C.FloatArray()
	uplo := linalg.ParamString(params.Uplo)
	trans := linalg.ParamString(params.Trans)
	//diag := linalg.ParamString(params.Diag)
	dsyr2k(uplo, trans, ind.N, ind.K, alpha, Aa[ind.OffsetA:], ind.LDa,
		Ba[ind.OffsetB:], ind.LDb, beta, Ca[ind.OffsetC:], ind.LDc)
	return
}

// See function Trmm.
func TrmmFloat(A, B *matrix.FloatMatrix, alpha float64, opts ...linalg.Option) (err error) {

	params, e := linalg.GetParameters(opts...)
	if e != nil {
		err = e
		return
	}
	ind := linalg.GetIndexOpts(opts...)
	err = check_level3_func(ind, ftrmm, A, B, nil, params)
	if err != nil {
		return
	}
	if ind.M == 0 || ind.N == 0 {
		return
	}
	Aa := A.FloatArray()
	Ba := B.FloatArray()
	uplo := linalg.ParamString(params.Uplo)
	transA := linalg.ParamString(params.TransA)
	side := linalg.ParamString(params.Side)
	diag := linalg.ParamString(params.Diag)
	dtrmm(side, uplo, transA, diag, ind.M, ind.N, alpha,
		Aa[ind.OffsetA:], ind.LDa, Ba[ind.OffsetB:], ind.LDb)
	return
}

// See function Trsm.
func TrsmFloat(A, B *matrix.FloatMatrix, alpha float64, opts ...linalg.Option) (err error) {

	params, e := linalg.GetParameters(opts...)
	if e != nil {
		err = e
		return
	}
	ind := linalg.GetIndexOpts(opts...)
	err = check_level3_func(ind, ftrsm, A, B, nil, params)
	if err != nil {
		return
	}
	if ind.N == 0 || ind.M == 0 {
		return
	}
	Aa := A.FloatArray()
	Ba := B.FloatArray()
	uplo := linalg.ParamString(params.Uplo)
	transA := linalg.ParamString(params.TransA)
	side := linalg.ParamString(params.Side)
	diag := linalg.ParamString(params.Diag)
	dtrsm(side, uplo, transA, diag, ind.M, ind.N, alpha,
		Aa[ind.OffsetA:], ind.LDa, Ba[ind.OffsetB:], ind.LDb)
	return
}

// Local Variables:
// tab-width: 4
// End:
