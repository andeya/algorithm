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

// Returns the Euclidean norm of a vector (returns ||x||_2).
//
// ARGUMENTS
//  X         float or complex matrix
//
// OPTIONS
//  n         integer.  If n<0, the default value of n is used.
//            The default value is equal to 1+(len(x)-offsetx-1)/incx or 0
//            if len(x) > offsetx+1
//  inc       positive integer
//  offset    nonnegative integer
//
func Nrm2(X matrix.Matrix, opts ...linalg.Option) (v matrix.Scalar) {
	v = matrix.FScalar(math.NaN())
	ind := linalg.GetIndexOpts(opts...)
	err := check_level1_func(ind, fnrm2, X, nil)
	if err != nil {
		return
	}
	if ind.Nx == 0 {
		return
	}
	switch X.(type) {
	case *matrix.ComplexMatrix:
		Xa := X.(*matrix.ComplexMatrix).ComplexArray()
		v = matrix.FScalar(dznrm2(ind.Nx, Xa[ind.OffsetX:], ind.IncX))
	case *matrix.FloatMatrix:
		Xa := X.(*matrix.FloatMatrix).FloatArray()
		v = matrix.FScalar(dnrm2(ind.Nx, Xa[ind.OffsetX:], ind.IncX))
	default:
		//err = onError("not implemented for parameter types", )
	}
	return
}

// Returns ||Re x||_1 + ||Im x||_1.
//
// ARGUMENTS
//  X       float or complex matrix
//
// OPTIONS
//  n       integer.  If n<0, the default value of n is used.
//          The default value is equal to n = 1+(len(x)-offset-1)/inc or 0 if
//          len(x) > offset+1
//  inc     positive integer
//  offset  nonnegative integer
//
func Asum(X matrix.Matrix, opts ...linalg.Option) (v matrix.Scalar) {
	v = matrix.FScalar(math.NaN())
	ind := linalg.GetIndexOpts(opts...)
	err := check_level1_func(ind, fasum, X, nil)
	if err != nil {
		return
	}
	if ind.Nx == 0 {
		return
	}
	switch X.(type) {
	case *matrix.ComplexMatrix:
		Xa := X.(*matrix.ComplexMatrix).ComplexArray()
		v = matrix.FScalar(dzasum(ind.Nx, Xa[ind.OffsetX:], ind.IncX))
	case *matrix.FloatMatrix:
		Xa := X.(*matrix.FloatMatrix).FloatArray()
		v = matrix.FScalar(dasum(ind.Nx, Xa[ind.OffsetX:], ind.IncX))
		//default:
		//	err = onError("not implemented for parameter types", )
	}
	return
}

// Returns Y = X^T*Y for real or complex X, Y.
//
// ARGUMENTS
//  X         float or complex matrix
//  Y         float or complex matrix.  Must have the same type as X.
//
// OPTIONS
//  n         integer.  If n<0, the default value of n is used.
//            The default value is equal to nx = 1+(len(x)-offsetx-1)/incx or 0 if
//            len(x) > offsetx+1.  If the default value is used, it must be equal to
//            ny = 1+(len(y)-offsetx-1)/|incy| or 0 if len(y) > offsety+1
//  incx      nonzero integer, [default=1]
//  incy      nonzero integer, [default=1]
//  offsetx   nonnegative integer, [default=0]
//  offsety   nonnegative integer, [default=0]
//
func Dotu(X, Y matrix.Matrix, opts ...linalg.Option) (v matrix.Scalar) {
	v = matrix.FScalar(math.NaN())
	//cv = cmplx.NaN()
	ind := linalg.GetIndexOpts(opts...)
	err := check_level1_func(ind, fdot, X, Y)
	if err != nil {
		return
	}
	if ind.Nx == 0 {
		return
	}
	sameType := matrix.EqualTypes(X, Y)
	if !sameType {
		err = onError("arrays not of same type")
		return
	}
	switch X.(type) {
	case *matrix.ComplexMatrix:
		Xa := X.(*matrix.ComplexMatrix).ComplexArray()
		Ya := Y.(*matrix.ComplexMatrix).ComplexArray()
		v = matrix.CScalar(zdotu(ind.Nx, Xa[ind.OffsetX:], ind.IncX, Ya[ind.OffsetY:], ind.IncY))
	case *matrix.FloatMatrix:
		Xa := X.(*matrix.FloatMatrix).FloatArray()
		Ya := Y.(*matrix.FloatMatrix).FloatArray()
		v = matrix.FScalar(ddot(ind.Nx, Xa[ind.OffsetX:], ind.IncX, Ya[ind.OffsetY:], ind.IncY))
		//default:
		//	err = onError("not implemented for parameter types", )
	}
	return
}

// Returns Y = X^H*Y for real or complex X, Y.
//
// ARGUMENTS
//  X         float or complex matrix
//  Y         float or complex matrix.  Must have the same type as X.
//
// OPTIONS
//  n         integer.  If n<0, the default value of n is used.
//            The default value is equal to nx = 1+(len(x)-offsetx-1)/incx or 0 if
//            len(x) > offsetx+1.  If the default value is used, it must be equal to
//            ny = 1+(len(y)-offsetx-1)/|incy| or 0 if len(y) > offsety+1
//  incx      nonzero integer [default=1]
//  incy      nonzero integer [default=1]
//  offsetx   nonnegative integer [default=0]
//  offsety   nonnegative integer [default=0]
//
func Dot(X, Y matrix.Matrix, opts ...linalg.Option) (v matrix.Scalar) {
	v = matrix.FScalar(math.NaN())
	//cv = cmplx.NaN()
	ind := linalg.GetIndexOpts(opts...)
	err := check_level1_func(ind, fdot, X, Y)
	if err != nil {
		return
	}
	if ind.Nx == 0 {
		return matrix.FScalar(0.0)
	}
	sameType := matrix.EqualTypes(X, Y)
	if !sameType {
		err = onError("arrays not of same type")
		return
	}
	switch X.(type) {
	case *matrix.ComplexMatrix:
		Xa := X.(*matrix.ComplexMatrix).ComplexArray()
		Ya := Y.(*matrix.ComplexMatrix).ComplexArray()
		v = matrix.CScalar(zdotc(ind.Nx, Xa[ind.OffsetX:], ind.IncX, Ya[ind.OffsetY:], ind.IncY))
	case *matrix.FloatMatrix:
		Xa := X.(*matrix.FloatMatrix).FloatArray()
		Ya := Y.(*matrix.FloatMatrix).FloatArray()
		v = matrix.FScalar(ddot(ind.Nx, Xa[ind.OffsetX:], ind.IncX, Ya[ind.OffsetY:], ind.IncY))
		//default:
		//	err = onError("not implemented for parameter types", )
	}
	return
}

// Interchanges two vectors (X <-> Y).
//
// ARGUMENTS
//  X         float or complex matrix
//  Y         float or complex matrix.  Must have the same type as X.
//
// OPTIONS
//  n         integer.  If n<0, the default value of n is used.
//            The default value is equal to 1+(len(x)-offsetx-1)/abs(incx) or
//            0 if len(x) > offsetx+1. Also if the default value is used,
//            it must be equal to 1+(len(y)-offsetx-1)/abs(incy) or 0 if
//            len(y) > offsety + 1.
//  incx      nonzero integer
//  incy      nonzero integer
//  offsetx   nonnegative integer
//  offsety   nonnegative integer;
//
func Swap(X, Y matrix.Matrix, opts ...linalg.Option) (err error) {
	ind := linalg.GetIndexOpts(opts...)
	err = check_level1_func(ind, fswap, X, Y)
	if err != nil {
		return
	}
	if ind.Nx == 0 {
		return
	}
	//if ind.Nx != ind.Ny {
	//	err = onError("arrays have unequal default lengths")
	//	return
	//}
	sameType := matrix.EqualTypes(X, Y)
	if !sameType {
		err = onError("arrays not same type")
		return
	}
	switch X.(type) {
	case *matrix.ComplexMatrix:
		Xa := X.(*matrix.ComplexMatrix).ComplexArray()
		Ya := Y.(*matrix.ComplexMatrix).ComplexArray()
		zswap(ind.Nx, Xa[ind.OffsetX:], ind.IncX, Ya[ind.OffsetY:], ind.IncY)
	case *matrix.FloatMatrix:
		Xa := X.(*matrix.FloatMatrix).FloatArray()
		Ya := Y.(*matrix.FloatMatrix).FloatArray()
		dswap(ind.Nx, Xa[ind.OffsetX:], ind.IncX, Ya[ind.OffsetY:], ind.IncY)
	default:
		err = onError("not implemented for parameter types")
	}
	return
}

// Copies a vector X to a vector Y (Y := X).
//
// ARGUMENTS
//  X         float or complex matrix
//  Y         float or complex matrix.  Must have the same type as X.
//
// OPTIONS
//  n         integer.  If n<0, the default value of n is used.
//            The default value is given by 1+(len(x)-offsetx-1)/incx or 0
//            if len(x) > offsetx+1
//  incx      nonzero integer
//  incy      nonzero integer
//  offsetx   nonnegative integer
//  offsety   nonnegative integer;
//
func Copy(X, Y matrix.Matrix, opts ...linalg.Option) (err error) {
	ind := linalg.GetIndexOpts(opts...)
	err = check_level1_func(ind, fcopy, X, Y)
	if err != nil {
		return
	}
	if ind.Nx == 0 {
		return
	}
	sameType := matrix.EqualTypes(X, Y)
	if !sameType {
		err = onError("arrays not same type")
		return
	}
	switch X.(type) {
	case *matrix.ComplexMatrix:
		Xa := X.(*matrix.ComplexMatrix).ComplexArray()
		Ya := Y.(*matrix.ComplexMatrix).ComplexArray()
		zcopy(ind.Nx, Xa[ind.OffsetX:], ind.IncX, Ya[ind.OffsetY:], ind.IncY)
	case *matrix.FloatMatrix:
		Xa := X.(*matrix.FloatMatrix).FloatArray()
		Ya := Y.(*matrix.FloatMatrix).FloatArray()
		dcopy(ind.Nx, Xa[ind.OffsetX:], ind.IncX, Ya[ind.OffsetY:], ind.IncY)
	default:
		err = onError("not implemented for parameter types")
	}
	return
}

// Scales a vector by a constant (X := alpha*X).
//
// ARGUMENTS
//  X         float or complex matrix
//  alpha     number (float or complex singleton matrix).  Complex alpha is only
//            allowed if X is complex.
//
// OPTIONS
//  n         integer.  If n<0, the default value of n is used.
//            The default value is equal to 1+(len(x)-offset-1)/inc or 0
//            if len(x) > offset+1.
//  inc       positive integer, default = 1
//  offset    nonnegative integer, default = 0
//
func Scal(X matrix.Matrix, alpha matrix.Scalar, opts ...linalg.Option) (err error) {
	ind := linalg.GetIndexOpts(opts...)
	err = check_level1_func(ind, fscal, X, nil)
	if err != nil {
		return
	}
	if ind.Nx == 0 {
		return
	}
	switch X.(type) {
	case *matrix.ComplexMatrix:
		Xa := X.(*matrix.ComplexMatrix).ComplexArray()
		cval := alpha.Complex()
		zscal(ind.Nx, cval, Xa[ind.OffsetX:], ind.IncX)
	case *matrix.FloatMatrix:
		Xa := X.(*matrix.FloatMatrix).FloatArray()
		rval := alpha.Float()
		if math.IsNaN(rval) {
			return onError("alpha not float value")
		}
		dscal(ind.Nx, rval, Xa[ind.OffsetX:], ind.IncX)
	default:
		err = onError("not implemented for parameter types")
	}
	return
}

// Calculate Y := alpha * X + Y. Y is set to new values.
// valid options: N, inc, incx, incy, offset, offsetx, offsety

// Constant times a vector plus a vector (Y := alpha*X+Y).
//
// ARGUMENTS
//   X         float or complex matrix
//   Y         float or complex matrix.  Must have the same type as X.
//   alpha     number (float or complex singleton matrix).  Complex alpha is only
//             allowed if x is complex.
//
// OPTIONS
//   n         integer.  If n<0, the default value of n is used.
//             The default value is equal to 1+(len(x)-offsetx-1)/incx
//             or 0 if  len(x) >= offsetx+1
//   incx      nonzero integer
//   incy      nonzero integer
//   offsetx   nonnegative integer
//   offsety   nonnegative integer;
//
func Axpy(X, Y matrix.Matrix, alpha matrix.Scalar, opts ...linalg.Option) (err error) {
	ind := linalg.GetIndexOpts(opts...)
	err = check_level1_func(ind, faxpy, X, Y)
	if err != nil {
		return
	}
	if ind.Nx == 0 {
		return
	}
	sameType := matrix.EqualTypes(X, Y)
	if !sameType {
		err = onError("arrays not same type")
		return
	}
	switch X.(type) {
	case *matrix.ComplexMatrix:
		Xa := X.(*matrix.ComplexMatrix).ComplexArray()
		Ya := Y.(*matrix.ComplexMatrix).ComplexArray()
		aval := alpha.Complex()
		if cmplx.IsNaN(aval) {
			return onError("alpha not complex value")
		}
		zaxpy(ind.Nx, aval, Xa[ind.OffsetX:],
			ind.IncX, Ya[ind.OffsetY:], ind.IncY)
	case *matrix.FloatMatrix:
		Xa := X.(*matrix.FloatMatrix).FloatArray()
		Ya := Y.(*matrix.FloatMatrix).FloatArray()
		aval := alpha.Float()
		if math.IsNaN(aval) {
			return onError("alpha not float value")
		}
		daxpy(ind.Nx, aval, Xa[ind.OffsetX:],
			ind.IncX, Ya[ind.OffsetY:], ind.IncY)
	default:
		err = onError("not implemented for parameter types")
	}
	return
}

// Local Variables:
// tab-width: 4
// End:
