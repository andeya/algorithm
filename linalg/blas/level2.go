// Copyright (c) Harri Rautila, 2012,2013

// This file is part of github.com/hrautila/linalg/blas package.
// It is free software, distributed under the terms of GNU Lesser General Public
// License Version 3, or any later version. See the COPYING tile included in this archive.

package blas

import (
	//"errors"
	"github.com/henrylee2cn/algorithm/linalg"
	"github.com/henrylee2cn/algorithm/matrix"
	"math"
	"math/cmplx"
)

/*
 General matrix-vector product.  (L2)

 Computes:
  Y := alpha*A*X + beta*Y,   if trans is NoTrans
  Y := alpha*A^T*X + beta*Y, if trans is Trans
  Y := beta*Y,               if n=0, m>0 and trans is NoTrans
  Y := beta*Y,               if n>0, m=0 and trans is Trans

 The matrix A is m by n. Returns immediately if n=0 and trans is PTrans,
 or if m=0 and trans is PNoTrans.

 ARGUMENTS
  A         float or complex m*n matrix
  X         float or complex n*1 matrix.
  Y         float or complex m*1 matrix.
  alpha     number (float or complex singleton matrix)
  beta      number (float or complex singleton matrix)

 OPTIONS
  trans     PNoTrans, PTrans, PConjTrans
  m         integer.  If negative, the default value is used.
  n         integer.  If negative, the default value is used.
  ldA       nonnegative integer.  ldA >= max(1,m). If zero, the default value is used.
  incx      nonzero integer
  incy      nonzero integer
  offsetA   nonnegative integer
  offsetx   nonnegative integer
  offsety   nonnegative integer
*/
func Gemv(A, X, Y matrix.Matrix, alpha, beta matrix.Scalar, opts ...linalg.Option) (err error) {

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
	if !matrix.EqualTypes(A, X, Y) {
		return onError("Parameters not of same type")
	}
	switch X.(type) {
	case *matrix.FloatMatrix:
		Xa := X.(*matrix.FloatMatrix).FloatArray()
		Ya := Y.(*matrix.FloatMatrix).FloatArray()
		Aa := A.(*matrix.FloatMatrix).FloatArray()
		aval := alpha.Float()
		bval := beta.Float()
		if math.IsNaN(aval) || math.IsNaN(bval) {
			return onError("alpha or beta not a number")
		}
		if params.Trans == linalg.PNoTrans && ind.N == 0 {
			dscal(ind.M, bval, Ya[ind.OffsetY:], ind.IncY)
		} else if params.Trans == linalg.PTrans && ind.M == 0 {
			dscal(ind.N, bval, Ya[ind.OffsetY:], ind.IncY)
		} else {
			trans := linalg.ParamString(params.Trans)
			dgemv(trans, ind.M, ind.N, aval, Aa[ind.OffsetA:],
				ind.LDa, Xa[ind.OffsetX:], ind.IncX, bval,
				Ya[ind.OffsetY:], ind.IncY)
		}
	case *matrix.ComplexMatrix:
		return onError("Not implemented yet for complx.Matrix")
	default:
		return onError("Unknown type, not implemented")
	}
	return
}

/*
 Matrix-vector product with a general banded matrix. (L2)

 Computes
   Y := alpha*A*X + beta*Y,   if trans = PNoTrans
   Y := alpha*A^T*X + beta*Y, if trans = PTrans
   Y := beta*y,               if n=0, m>0, and trans = PNoTrans
   Y := beta*y,               if n>0, m=0, and trans = PTrans

 The matrix A is m by n with upper bandwidth ku and lower bandwidth kl.
 Returns immediately if n=0 and trans is 'Trans', or if m=0 and trans is 'N'.


 ARGUMENTS
   X         float n*1 matrix.
   Y         float m*1 matrix
   A         float m*n matrix.
   alpha     number (float).
   beta      number (float).

 OPTIONS
   trans     NoTrans or Trans
   m         nonnegative integer, default A.Rows()
   kl        nonnegative integer
   n         nonnegative integer.  If negative, the default value is used.
   ku        nonnegative integer.  If negative, the default value is used.
   ldA       positive integer.  ldA >= kl+ku+1. If zero, the default value is used.
   incx      nonzero integer, default =1
   incy      nonzero integer, default =1
   offsetA   nonnegative integer, default =0
   offsetx   nonnegative integer, default =0
   offsety   nonnegative integer, default =0

*/
func Gbmv(A, X, Y matrix.Matrix, alpha, beta matrix.Scalar, opts ...linalg.Option) (err error) {

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
	if !matrix.EqualTypes(A, X, Y) {
		return onError("Parameters not of same type")
	}
	switch X.(type) {
	case *matrix.FloatMatrix:
		Xa := X.(*matrix.FloatMatrix).FloatArray()
		Ya := Y.(*matrix.FloatMatrix).FloatArray()
		Aa := A.(*matrix.FloatMatrix).FloatArray()
		aval := alpha.Float()
		bval := beta.Float()
		if math.IsNaN(aval) || math.IsNaN(bval) {
			return onError("alpha or beta not a number")
		}
		if params.Trans == linalg.PNoTrans && ind.N == 0 {
			dscal(ind.M, bval, Ya[ind.OffsetY:], ind.IncY)
		} else if params.Trans == linalg.PTrans && ind.M == 0 {
			dscal(ind.N, bval, Ya[ind.OffsetY:], ind.IncY)
		} else {
			trans := linalg.ParamString(params.Trans)
			dgbmv(trans, ind.M, ind.N, ind.Kl, ind.Ku,
				aval, Aa[ind.OffsetA:], ind.LDa, Xa[ind.OffsetX:], ind.IncX,
				bval, Ya[ind.OffsetY:], ind.IncY)
		}
	case *matrix.ComplexMatrix:
		return onError("Not implemented yet for complx.Matrix")
	default:
		return onError("Unknown type, not implemented")
	}
	return
}

/*
 Matrix-vector product with a real symmetric matrix. (L2)

 Computes with A real symmetric of order n.
  Y := alpha*A*X + beta*Y

 ARGUMENTS
  A         float or complex n*n matrix
  X         float or complex n*1 matrix
  Y         float or complex n*1 matrix
  alpha     number (float or complex singleton matrix)
  beta      number (float or complex singleton matrix)

 OPTIONS
  uplo      PLower or PUpper
  n         integer.  If negative, the default value is used.
            If the default value is used, we require that A.Rows()=A.Cols().
  ldA       nonnegative integer.  ldA >= max(1,n). If zero, the default value is used.
  incx      nonzero integer
  incy      nonzero integer
  offsetA   nonnegative integer
  offsetx   nonnegative integer
  offsety   nonnegative integer
*/
func Symv(A, X, Y matrix.Matrix, alpha, beta matrix.Scalar, opts ...linalg.Option) (err error) {

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
	if !matrix.EqualTypes(A, X, Y) {
		return onError("Parameters not of same type")
	}
	switch X.(type) {
	case *matrix.FloatMatrix:
		Xa := X.(*matrix.FloatMatrix).FloatArray()
		Ya := Y.(*matrix.FloatMatrix).FloatArray()
		Aa := A.(*matrix.FloatMatrix).FloatArray()
		aval := alpha.Float()
		bval := beta.Float()
		if math.IsNaN(aval) || math.IsNaN(bval) {
			return onError("alpha or beta not a number")
		}
		uplo := linalg.ParamString(params.Uplo)
		dsymv(uplo, ind.N, aval, Aa[ind.OffsetA:], ind.LDa,
			Xa[ind.OffsetX:], ind.IncX,
			bval, Ya[ind.OffsetY:], ind.IncY)
	case *matrix.ComplexMatrix:
		return onError("Symv not possible for ComplexMatrix")
	default:
		return onError("Unknown type, not implemented")
	}
	return
}

/*
 Matrix-vector product with a real symmetric or complex Hermitian matrix. (L2)

 Hemv(A, X, Y, alpha=1.0, beta=0.0, uplo=PLower, n=A.Rows,
 ldA=max(1,A.Rows), incx=1, incy=1, offsetA=0, offsetx=0, offsety=0)

 Computes
  Y := alpha*A*X + beta*Y  with A real symmetric of order n.

 ARGUMENTS
  A         float or complex n*n matrix
  X         float or complex n*1 matrix
  Y         float or complex n*1 matrix
  alpha     number (float or complex singleton matrix)
  beta      number (float or complex singleton matrix)

 OPTIONS
  uplo      PLower or PUpper
  n         integer.  If negative, the default value is used.
            If the default value is used, we require that  A.Rows()=A.Cols().
  ldA       nonnegative integer.  ldA >= max(1,n). If zero, the default value is used.
  incx      nonzero integer
  incy      nonzero integer
  offsetA   nonnegative integer
  offsetx   nonnegative integer
  offsety   nonnegative integer
*/
func Hemv(A, X, Y matrix.Matrix, alpha, beta matrix.Scalar, opts ...linalg.Option) (err error) {

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
	if !matrix.EqualTypes(A, X, Y) {
		return onError("Parameters not of same type")
	}
	switch X.(type) {
	case *matrix.FloatMatrix:
		Xa := X.(*matrix.FloatMatrix).FloatArray()
		Ya := Y.(*matrix.FloatMatrix).FloatArray()
		Aa := A.(*matrix.FloatMatrix).FloatArray()
		aval := alpha.Float()
		bval := beta.Float()
		if math.IsNaN(aval) || math.IsNaN(bval) {
			return onError("alpha or beta not a number")
		}
		uplo := linalg.ParamString(params.Uplo)
		dsymv(uplo, ind.N, aval, Aa[ind.OffsetA:], ind.LDa,
			Xa[ind.OffsetX:], ind.IncX,
			bval, Ya[ind.OffsetY:], ind.IncY)

	case *matrix.ComplexMatrix:
		Xa := X.(*matrix.ComplexMatrix).ComplexArray()
		Ya := Y.(*matrix.ComplexMatrix).ComplexArray()
		Aa := A.(*matrix.ComplexMatrix).ComplexArray()
		aval := alpha.Complex()
		/*
			if alpha != nil {
				aval = alpha.ComplexValue()
				if cmplx.IsNaN(aval) {
					if ! math.IsNaN(alpha.FloatValue()) {
						aval = complex(alpha.FloatValue(), 0)
					} else {
						return onError("alpha not a number")
					}
				}
			}
		*/
		bval := beta.Complex()
		/*
			if beta != nil {
				bval := beta.ComplexValue()
				if cmplx.IsNaN(bval) {
					if ! math.IsNaN(beta.FloatValue()) {
						bval = complex(beta.FloatValue(), 0)
					} else {
						return onError("beta not a number")
					}
				}
			}
		*/
		uplo := linalg.ParamString(params.Uplo)
		zhemv(uplo, ind.N, aval, Aa[ind.OffsetA:], ind.LDa,
			Xa[ind.OffsetX:], ind.IncX,
			bval, Ya[ind.OffsetY:], ind.IncY)
	default:
		return onError("Unknown type, not implemented")
	}
	return
}

/*
 Matrix-vector product with a real symmetric band matrix. (L2)

 Sbmv(A, X, Y, alpha=1.0, beta=0.0, uplo=PLower, n=A.Cols,
 k=-1, ldA=A.Rows, incx=1, incy=1, offsetA=0, offsetx=0, offsety=0)

 Computes with A real symmetric and  banded of order n and with bandwidth k.
   Y := alpha*A*X + beta*Y

 ARGUMENTS
   A         float or complex n*n matrix
   X         float or complex n*1 matrix
   Y         float or complex n*1 matrix
   alpha     number (float or complex singleton matrix)
   beta      number (float or complex singleton matrix)

 OPTIONS
   uplo      PLower or PUpper
   n         integer.  If negative, the default value is used.
   k         integer.  If negative, the default value is used. The default value is
 		     k = max(0,A.Rows()-1).
   ldA       nonnegative integer.  ldA >= k+1. If zero, the default vaule is used.
   incx      nonzero integer
   incy      nonzero integer
   offsetA   nonnegative integer
   offsetx   nonnegative integer
   offsety   nonnegative integer

*/
func Sbmv(A, X, Y matrix.Matrix, alpha, beta matrix.Scalar, opts ...linalg.Option) (err error) {

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
	if !matrix.EqualTypes(A, X, Y) {
		return onError("Parameters not of same type")
	}
	switch X.(type) {
	case *matrix.FloatMatrix:
		Xa := X.(*matrix.FloatMatrix).FloatArray()
		Ya := Y.(*matrix.FloatMatrix).FloatArray()
		Aa := A.(*matrix.FloatMatrix).FloatArray()
		aval := alpha.Float()
		bval := beta.Float()
		if math.IsNaN(aval) || math.IsNaN(bval) {
			return onError("alpha or beta not a number")
		}
		uplo := linalg.ParamString(params.Uplo)
		dsbmv(uplo, ind.N, ind.K, aval, Aa[ind.OffsetA:], ind.LDa,
			Xa[ind.OffsetX:], ind.IncX, bval, Ya[ind.OffsetY:], ind.IncY)

	case *matrix.ComplexMatrix:
		return onError("Sbmv not possible for complx.Matrix")
	default:
		return onError("Unknown type, not implemented")
	}
	return
}

/*
 Matrix-vector product with a real symmetric or complex hermitian band matrix.

 Computes with A real symmetric and  banded of order n and with bandwidth k.
  Y := alpha*A*X + beta*Y

 ARGUMENTS
  A         float or complex n*n matrix
  X         float or complex n*1 matrix
  Y         float or complex n*1 matrix
  alpha     number (float or complex singleton matrix)
  beta      number (float or complex singleton matrix)

 OPTIONS
  uplo      PLower or PUpper
  n         integer.  If negative, the default value is used.
  k         integer.  If negative, the default value is used.
            The default value is k = max(0,A.Rows()-1).
  ldA       nonnegative integer.  ldA >= k+1.
            If zero, the default vaule is used.
  incx      nonzero integer
  incy      nonzero integer
  offsetA   nonnegative integer
  offsetx   nonnegative integer
  offsety   nonnegative integer

*/
func Hbmv(A, X, Y matrix.Matrix, alpha, beta matrix.Scalar, opts ...linalg.Option) (err error) {

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
	if !matrix.EqualTypes(A, X, Y) {
		return onError("Parameters not of same type")
	}
	switch X.(type) {
	case *matrix.FloatMatrix:
		Xa := X.(*matrix.FloatMatrix).FloatArray()
		Ya := Y.(*matrix.FloatMatrix).FloatArray()
		Aa := A.(*matrix.FloatMatrix).FloatArray()
		aval := alpha.Float()
		bval := beta.Float()
		if math.IsNaN(aval) || math.IsNaN(bval) {
			return onError("alpha or beta not a number")
		}
		uplo := linalg.ParamString(params.Uplo)
		dsbmv(uplo, ind.N, ind.K, aval, Aa[ind.OffsetA:], ind.LDa,
			Xa[ind.OffsetX:], ind.IncX, bval, Ya[ind.OffsetY:], ind.IncY)

	case *matrix.ComplexMatrix:
		Xa := X.(*matrix.ComplexMatrix).ComplexArray()
		Ya := Y.(*matrix.ComplexMatrix).ComplexArray()
		Aa := A.(*matrix.ComplexMatrix).ComplexArray()
		aval := alpha.Complex()
		bval := beta.Complex()
		uplo := linalg.ParamString(params.Uplo)
		zhbmv(uplo, ind.N, ind.K, aval, Aa[ind.OffsetA:], ind.LDa,
			Xa[ind.OffsetX:], ind.IncX, bval, Ya[ind.OffsetY:], ind.IncY)
		//zhbmv(uplo, ind.N, aval, Aa[ind.OffsetA:], ind.LDa,
		//	Xa[ind.OffsetX:], ind.IncX,
		//	bval, Ya[ind.OffsetY:], ind.IncY)
	default:
		return onError("Unknown type, not implemented")
	}
	return
}

/*
 Matrix-vector product with a triangular matrix. (L2)

 Trmv(A, x, uplo=PLower, trans=PNoTrans, diag=PNonUnit, n=A.Rows,
 ldA=max(1,A.Rows), incx=1, offsetA=0, offsetx=0)

 COMPUTES
  X := A*X,   if trans is PNoTrans
  X := A^T*X, if trans is PTrans
  X := A^H*X, if trans is PConjTrans

 A is triangular of order n.

 ARGUMENTS
  A         float or complex matrix
  X         float or complex  matrix.  Must have the same type as A.

 OPTIONS
  uplo      PLower or PUpper
  trans     PNoTrans or PTrans
  diag      PNonUnit or PUnit
  n         integer.  If negative, the default value is used.
            If the default value is used, we require that A.Rows = A.Cols.
  ldA       nonnegative integer.  ldA >= max(1,n).
            If zero the default value is used.
  incx      nonzero integer, default=1
  offsetA   nonnegative integer, default=0
  offsetx   nonnegative integer, default=0

*/
func Trmv(A, X matrix.Matrix, opts ...linalg.Option) (err error) {

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
	if !matrix.EqualTypes(A, X) {
		return onError("Parameters not of same type")
	}
	switch X.(type) {
	case *matrix.FloatMatrix:
		Xa := X.(*matrix.FloatMatrix).FloatArray()
		Aa := A.(*matrix.FloatMatrix).FloatArray()
		uplo := linalg.ParamString(params.Uplo)
		trans := linalg.ParamString(params.Trans)
		diag := linalg.ParamString(params.Diag)
		dtrmv(uplo, trans, diag, ind.N,
			Aa[ind.OffsetA:], ind.LDa, Xa[ind.OffsetX:], ind.IncX)
	case *matrix.ComplexMatrix:
		return onError("Not implemented yet for complx.Matrix")
	default:
		return onError("Unknown type, not implemented")
	}
	return
}

/*
 Matrix-vector product with a triangular band matrix. (L2)

 Tbmv(A, x, uplo=PLower, trans=PNoTrans, diag=PNonUnit, n=A.Cols,
 k=max(0,A.Rows-1), ldA=A.size[0], incx=1, offsetA=0, offsetx=0)

 COMPUTES
  X := A*X,   if trans is PNoTrans
  X := A^T*X, if trans is PTrans
  X := A^H*X, if trans is PConjTrans

 A is banded triangular of order n and with bandwith k.

 ARGUMENTS
  A         float or complex matrix
  X         float or complex  matrix.  Must have the same type as A.

 OPTIONS
  uplo      PLower or PUpper
  trans     PNoTrans or PTrans
  diag      PNonUnit or PUnit
  n         nonnegative integer.  If negative, the default value is used.
  k         nonnegative integer.  If negative, the default value is used.
  ldA       nonnegative integer.  lda >= 1+k.
            If zero the default value is used.
  incx      nonzero integer
  offsetA   nonnegative integer
  offsetx   nonnegative integer

*/
func Tbmv(A, X matrix.Matrix, opts ...linalg.Option) (err error) {

	var params *linalg.Parameters
	if !matrix.EqualTypes(A, X) {
		err = onError("Parameters not of same type")
		return
	}
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
	switch X.(type) {
	case *matrix.FloatMatrix:
		Xa := X.(*matrix.FloatMatrix).FloatArray()
		Aa := A.(*matrix.FloatMatrix).FloatArray()
		uplo := linalg.ParamString(params.Uplo)
		trans := linalg.ParamString(params.Trans)
		diag := linalg.ParamString(params.Diag)
		dtbmv(uplo, trans, diag, ind.N, ind.K,
			Aa[ind.OffsetA:], ind.LDa, Xa[ind.OffsetX:], ind.IncX)
	case *matrix.ComplexMatrix:
		return onError("Not implemented yet for ComplexMatrix")
	default:
		return onError("Unknown type, not implemented")
	}
	return
}

/*
 Solution of a triangular set of equations with one righthand side. (L2)

 Trsv(A, x, uplo=PLower, trans=PNoTrans, diag=PNonUnit, n=A.Cols,
 ldA=max(1,A.Rows), incx=1, offsetA=0, offsetx=0)

 PURPOSE
  X := A^{-1}*X, if trans is PNoTrans
  X := A^{-T}*X, if trans is PTrans
  X := A^{-H}*X, if trans is PConjTrans

 A is triangular of order n.  The code does not verify whether A is nonsingular.

 ARGUMENTS
  A         float or complex m*n matrix.
  X         float or complex n*1 matrix. Must have the same type as A.

 OPTIONS
  uplo      PLower   or PUpper
  trans     PNoTrans or PTrans
  diag      PNoNUnit or PUnit
  n         integer.  If negative, the default value is used.
            If the default value is used, we require that A.rows = A.cols.
  ldA       nonnegative integer.  ldA >= max(1,n). If zero, the default value is used.
  incx      nonzero integer
  offsetA   nonnegative integer
  offsetx   nonnegative integer
*/
func Trsv(A, X matrix.Matrix, opts ...linalg.Option) (err error) {

	var params *linalg.Parameters
	if !matrix.EqualTypes(A, X) {
		err = onError("Parameters not of same type")
		return
	}
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
	switch X.(type) {
	case *matrix.FloatMatrix:
		Xa := X.(*matrix.FloatMatrix).FloatArray()
		Aa := A.(*matrix.FloatMatrix).FloatArray()
		uplo := linalg.ParamString(params.Uplo)
		trans := linalg.ParamString(params.Trans)
		diag := linalg.ParamString(params.Diag)
		dtrsv(uplo, trans, diag, ind.N,
			Aa[ind.OffsetA:], ind.LDa, Xa[ind.OffsetX:], ind.IncX)
	case *matrix.ComplexMatrix:
		return onError("Not implemented yet for complx.Matrix")
	default:
		return onError("Unknown type, not implemented")
	}
	return
}

/*
 Solution of a triangular and banded set of equations.

 Tbsv(A, X, uplo=PLower, trans=PNoTrans, diag=PNonDiag, n=A.Cols,
 k=max(0,A.Rows-1), ldA=A.size[0], incx=1, offsetA=0, offsetx=0)

PURPOSE
  X := A^{-1}*X, if trans is PNoTrans
  X := A^{-T}*X, if trans is PTrans
  X := A^{-H}*X, if trans is PConjTrans

 A is banded triangular of order n and with bandwidth k.

 ARGUMENTS
  A         float or complex m*k matrix.
  X         float or complex k*1 matrix. Must have the same type as A.

 OPTIONS
  uplo      PLower   or PUpper
  trans     PNoTrans, PTrans or PConjTrans
  diag      PNoNUnit or PUnit
  n         nonnegative integer.  If negative, the default value is used.
  k         nonnegative integer.  If negative, the default value is used.
  ldA       nonnegative integer.  ldA >= 1+k.
            If zero the default value is used.
  incx      nonzero integer
  offsetA   nonnegative integer
  offsetx   nonnegative integer;
*/
func Tbsv(A, X matrix.Matrix, opts ...linalg.Option) (err error) {

	var params *linalg.Parameters
	if !matrix.EqualTypes(A, X) {
		err = onError("Parameters not of same type")
		return
	}
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
	switch X.(type) {
	case *matrix.FloatMatrix:
		Xa := X.(*matrix.FloatMatrix).FloatArray()
		Aa := A.(*matrix.FloatMatrix).FloatArray()
		uplo := linalg.ParamString(params.Uplo)
		trans := linalg.ParamString(params.Trans)
		diag := linalg.ParamString(params.Diag)
		dtbsv(uplo, trans, diag, ind.N, ind.K,
			Aa[ind.OffsetA:], ind.LDa, Xa[ind.OffsetX:], ind.IncX)
	case *matrix.ComplexMatrix:
		return onError("Not implemented yet for complx.Matrix")
	default:
		return onError("Unknown type, not implemented")
	}
	return
}

/*
 General rank-1 update. (L2)

 Ger(X, Y, A, alpha=1.0, m=A.Rows, n=A.Cols, incx=1,
 incy=1, ldA=max(1,A.Rows), offsetx=0, offsety=0, offsetA=0)

 COMPUTES
  A := A + alpha*X*Y^H with A m*n, real or complex.

 ARGUMENTS
  X         float or complex matrix.
  Y         float or complex matrix. Must have the same type as X.
  A         float or complex matrix. Must have the same type as X.
  alpha     number (float or complex singleton matrix).

 OPTIONS
  m         integer.  If negative, the default value is used.
  n         integer.  If negative, the default value is used.
  incx      nonzero integer
  incy      nonzero integer
  ldA       nonnegative integer.  ldA >= max(1,m).
            If zero, the default value is used.
  offsetx   nonnegative integer
  offsety   nonnegative integer
  offsetA   nonnegative integer;

*/
func Ger(X, Y, A matrix.Matrix, alpha matrix.Scalar, opts ...linalg.Option) (err error) {

	var params *linalg.Parameters
	if !matrix.EqualTypes(A, X, Y) {
		err = onError("Parameters not of same type")
		return
	}
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
	switch X.(type) {
	case *matrix.FloatMatrix:
		Xa := X.(*matrix.FloatMatrix).FloatArray()
		Ya := Y.(*matrix.FloatMatrix).FloatArray()
		Aa := A.(*matrix.FloatMatrix).FloatArray()
		aval := alpha.Float()
		if math.IsNaN(aval) {
			return onError("alpha not a number")
		}
		dger(ind.M, ind.N, aval, Xa[ind.OffsetX:], ind.IncX,
			Ya[ind.OffsetY:], ind.IncY, Aa[ind.OffsetA:], ind.LDa)

	case *matrix.ComplexMatrix:
		Xa := X.(*matrix.ComplexMatrix).ComplexArray()
		Ya := Y.(*matrix.ComplexMatrix).ComplexArray()
		Aa := A.(*matrix.ComplexMatrix).ComplexArray()
		aval := alpha.Complex()
		if cmplx.IsNaN(aval) {
			return onError("alpha not a number")
		}
		zgerc(ind.M, ind.N, aval, Xa[ind.OffsetX:], ind.IncX,
			Ya[ind.OffsetY:], ind.IncY, Aa[ind.OffsetA:], ind.LDa)

	default:
		return onError("Unknown type, not implemented")
	}
	return
}

/*
 General rank-1 update. (L2)

 Geru(X, Y, A, alpha=1.0, m=A.Rows, n=A.Cols, incx=1,
 incy=1, ldA=max(1,A.Rows), offsetx=0, offsety=0, offsetA=0)

 COMPUTES
  A := A + alpha*X*Y^T with A m*n, real or complex.

 ARGUMENTS
  X         float or complex matrix.
  Y         float or complex matrix. Must have the same type as X.
  A         float or complex matrix. Must have the same type as X.
  alpha     number (float or complex singleton matrix).

 OPTIONS
  m         integer.  If negative, the default value is used.
  n         integer.  If negative, the default value is used.
  incx      nonzero integer
  incy      nonzero integer
  ldA       nonnegative integer.  ldA >= max(1,m).
            If zero, the default value is used.
  offsetx   nonnegative integer
  offsety   nonnegative integer
  offsetA   nonnegative integer;

*/
func Geru(X, Y, A matrix.Matrix, alpha matrix.Scalar, opts ...linalg.Option) (err error) {

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
	if ind.M == 0 || ind.N == 0 {
		return
	}
	if !matrix.EqualTypes(A, X, Y) {
		return onError("Parameters not of same type")
	}
	switch X.(type) {
	case *matrix.FloatMatrix:
		Xa := X.(*matrix.FloatMatrix).FloatArray()
		Ya := Y.(*matrix.FloatMatrix).FloatArray()
		Aa := A.(*matrix.FloatMatrix).FloatArray()
		aval := alpha.Float()
		if math.IsNaN(aval) {
			return onError("alpha not a number")
		}
		dger(ind.M, ind.N, aval, Xa[ind.OffsetX:], ind.IncX,
			Ya[ind.OffsetY:], ind.IncY, Aa[ind.OffsetA:], ind.LDa)

	case *matrix.ComplexMatrix:
		Xa := X.(*matrix.ComplexMatrix).ComplexArray()
		Ya := Y.(*matrix.ComplexMatrix).ComplexArray()
		Aa := A.(*matrix.ComplexMatrix).ComplexArray()
		aval := alpha.Complex()
		if cmplx.IsNaN(aval) {
			return onError("alpha not a number")
		}
		zgeru(ind.M, ind.N, aval, Xa[ind.OffsetX:], ind.IncX,
			Ya[ind.OffsetY:], ind.IncY, Aa[ind.OffsetA:], ind.LDa)

	default:
		return onError("Unknown type, not implemented")
	}
	return
}

/*
 Symmetric rank-1 update. (L2)

 syr(X, A, alpha=1.0, uplo=PLower, n=A.Rows, incx=1,
 ldA=max(1,A.Rows), offsetx=0, offsetA=0)

 COMPUTES
  A := A + alpha*X*X^T

 A real symmetric matrix of order n.

 ARGUMENTS
  X         float or complex matrix.
  A         float or complex matrix.
  alpha     real number

 OPTIONS:
  uplo      PLower or PUpper
  n         integer.  If negative, the default value is used.
  incx      nonzero integer
  ldA       nonnegative integer.  ldA >= max(1,n).
            If zero, the default value is used.
  offsetx   nonnegative integer
  offsetA   nonnegative integer;
*/
func Syr(X, A matrix.Matrix, alpha matrix.Scalar, opts ...linalg.Option) (err error) {

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
	if !matrix.EqualTypes(A, X) {
		return onError("Parameters not of same type")
	}
	switch X.(type) {
	case *matrix.FloatMatrix:
		Xa := X.(*matrix.FloatMatrix).FloatArray()
		Aa := A.(*matrix.FloatMatrix).FloatArray()
		uplo := linalg.ParamString(params.Uplo)
		aval := alpha.Float()
		if math.IsNaN(aval) {
			return onError("alpha not a number")
		}
		dsyr(uplo, ind.N, aval, Xa[ind.OffsetX:], ind.IncX,
			Aa[ind.OffsetA:], ind.LDa)
	case *matrix.ComplexMatrix:
		return onError("Syr not possible for complx.Matrix")
	default:
		return onError("Unknown type, not implemented")
	}
	return
}

/*
 Symmetric rank-1 update. (L2)

 Her(X, A, alpha, uplo=PLowoer, n=A.Rows, incx=1,
 ldA=max(1,A.Rows), offsetx=0, offsetA=0)

 COMPUTES
  A := A + alpha*X*X^H

 A real symmetric or complex hermitian matrix of order n.

 ARGUMENTS
  X         float or complex matrix.
  A         float or complex matrix.
  alpha     number (float or complex singleton matrix)

 OPTIONS:
  uplo      PLower or PUpper
  n         integer.  If negative, the default value is used.
  incx      nonzero integer
  ldA       nonnegative integer.  ldA >= max(1,n).
            If zero, the default value is used.
  offsetx   nonnegative integer
  offsetA   nonnegative integer;

*/
func Her(X, A matrix.Matrix, alpha matrix.Scalar, opts ...linalg.Option) (err error) {

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
	if !matrix.EqualTypes(A, X) {
		return onError("Parameters not of same type")
	}
	switch X.(type) {
	case *matrix.FloatMatrix:
		Xa := X.(*matrix.FloatMatrix).FloatArray()
		Aa := A.(*matrix.FloatMatrix).FloatArray()
		uplo := linalg.ParamString(params.Uplo)
		aval := alpha.Float()
		if math.IsNaN(aval) {
			return onError("alpha not a number")
		}
		dsyr(uplo, ind.N, aval, Xa[ind.OffsetX:], ind.IncX,
			Aa[ind.OffsetA:], ind.LDa)
	case *matrix.ComplexMatrix:
		Xa := X.(*matrix.ComplexMatrix).ComplexArray()
		Aa := A.(*matrix.ComplexMatrix).ComplexArray()
		uplo := linalg.ParamString(params.Uplo)
		aval := alpha.Float()
		if math.IsNaN(aval) {
			return onError("alpha not a number")
		}
		zher(uplo, ind.N, aval, Xa[ind.OffsetX:], ind.IncX,
			Aa[ind.OffsetA:], ind.LDa)
	default:
		return onError("Unknown type, not implemented")
	}
	return
}

/*
 Symmetric rank-2 update.
 syr2(x, y, A, uplo='L', alpha=1.0, n=A.size[0], incx=1, incy=1,
     ldA=max(1,A.size[0]), offsetx=0, offsety=0, offsetA=0)
 PURPOSE
 Computes A := A + alpha*(x*y^T + y*x^T) with A real symmetric matrix of order n.
 ARGUMENTS
 x         float matrix
 y         float matrix
 A         float matrix
 alpha     real number (int or float)

 OPTIONS
 uplo      'L' or 'U'
 n         integer.  If negative, the default value is used.
 incx      nonzero integer
 incy      nonzero integer
 ldA       nonnegative integer.  ldA >= max(1,n).
           If zero the default value is used.
 offsetx   nonnegative integer
 offsety   nonnegative integer
 offsetA   nonnegative integer;
*/
func Syr2(X, Y, A matrix.Matrix, alpha matrix.Scalar, opts ...linalg.Option) (err error) {

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
	if !matrix.EqualTypes(A, X, Y) {
		return onError("Parameters not of same type")
	}
	switch X.(type) {
	case *matrix.FloatMatrix:
		Xa := X.(*matrix.FloatMatrix).FloatArray()
		Ya := X.(*matrix.FloatMatrix).FloatArray()
		Aa := A.(*matrix.FloatMatrix).FloatArray()
		aval := alpha.Float()
		if math.IsNaN(aval) {
			return onError("alpha not a number")
		}
		uplo := linalg.ParamString(params.Uplo)
		dsyr2(uplo, ind.N, aval, Xa[ind.OffsetX:], ind.IncX,
			Ya[ind.OffsetY:], ind.IncY,
			Aa[ind.OffsetA:], ind.LDa)
	case *matrix.ComplexMatrix:
		return onError("Not implemented yet for complx.Matrix")
	default:
		return onError("Unknown type, not implemented")
	}
	return
}

/*
 Symmetric rank-2 update.
 her2(x, y, A, uplo='L', alpha=1.0, n=A.size[0], incx=1, incy=1,
     ldA=max(1,A.size[0]), offsetx=0, offsety=0, offsetA=0)
 PURPOSE
 Computes A := A + alpha*(x*y^T + y*x^T) with A real symmetric or complex hermitian
 matix of order n.

 ARGUMENTS
 x         float or complex matrix
 y         float or complex matrix
 A         float or complex matrix
 alpha     float or complex singleton value

 OPTIONS:
 uplo      'L' or 'U'
 n         integer.  If negative, the default value is used.
 incx      nonzero integer
 incy      nonzero integer
 ldA       nonnegative integer.  ldA >= max(1,n).
           If zero the default value is used.
 offsetx   nonnegative integer
 offsety   nonnegative integer
 offsetA   nonnegative integer;
*/
func Her2(X, Y, A matrix.Matrix, alpha matrix.Scalar, opts ...linalg.Option) (err error) {

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
	if !matrix.EqualTypes(A, X, Y) {
		return onError("Parameters not of same type")
	}
	switch X.(type) {
	case *matrix.FloatMatrix:
		Xa := X.(*matrix.FloatMatrix).FloatArray()
		Ya := X.(*matrix.FloatMatrix).FloatArray()
		Aa := A.(*matrix.FloatMatrix).FloatArray()
		aval := alpha.Float()
		if math.IsNaN(aval) {
			return onError("alpha not a number")
		}
		uplo := linalg.ParamString(params.Uplo)
		dsyr2(uplo, ind.N, aval, Xa[ind.OffsetX:], ind.IncX,
			Ya[ind.OffsetY:], ind.IncY,
			Aa[ind.OffsetA:], ind.LDa)
	case *matrix.ComplexMatrix:
		Xa := X.(*matrix.ComplexMatrix).ComplexArray()
		Ya := X.(*matrix.ComplexMatrix).ComplexArray()
		Aa := A.(*matrix.ComplexMatrix).ComplexArray()
		aval := alpha.Complex()
		if cmplx.IsNaN(aval) {
			return onError("alpha not a number")
		}
		uplo := linalg.ParamString(params.Uplo)
		zher2(uplo, ind.N, aval, Xa[ind.OffsetX:], ind.IncX,
			Ya[ind.OffsetY:], ind.IncY,
			Aa[ind.OffsetA:], ind.LDa)
		//zher(uplo, ind.N, aval, Xa[ind.OffsetX:], ind.IncX,
		//	Aa[ind.OffsetA:], ind.LDa)
	default:
		return onError("Unknown type, not implemented")
	}
	return
}

// Local Variables:
// tab-width: 4
// End:
