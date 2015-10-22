// Copyright (c) Harri Rautila, 2012,2013

// This file is part of github.com/hrautila/linalg/blas package.
// It is free software, distributed under the terms of GNU Lesser General Public
// License Version 3, or any later version. See the COPYING tile included in this archive.

package blas

import (
	"github.com/henrylee2cn/algorithm/linalg"
	"github.com/henrylee2cn/algorithm/matrix"
	//"errors"
)

// Calculate norm2(X) for ComplexMatrix. Valid options are: offset, inc, N.

// See function Nrm2.
func Nrm2Complex(X *matrix.ComplexMatrix, opts ...linalg.Option) (v float64, err error) {
	v = 0.0
	ind := linalg.GetIndexOpts(opts...)
	err = check_level1_func(ind, fnrm2, X, nil)
	if err != nil {
		return
	}
	if ind.Nx == 0 {
		return
	}
	Xa := X.ComplexArray()
	v = dznrm2(ind.Nx, Xa[ind.OffsetX:], ind.IncX)
	return
}

// Calculate sum(X) for complex matrix. Valid options are: offset, inc, N.

// See function Asum.
func AsumComplex(X *matrix.ComplexMatrix, opts ...linalg.Option) (v float64, err error) {
	v = 0.0
	ind := linalg.GetIndexOpts(opts...)
	err = check_level1_func(ind, fasum, X, nil)
	if err != nil {
		return
	}
	if ind.Nx == 0 {
		return
	}
	Xa := X.ComplexArray()
	v = dzasum(ind.Nx, Xa[ind.OffsetX:], ind.IncX)
	return
}

// Calculate X.T * Y for complex matrix.

// See function Dot.
func DotuComplex(X, Y *matrix.ComplexMatrix, opts ...linalg.Option) (v complex128, err error) {
	v = 0.0
	ind := linalg.GetIndexOpts(opts...)
	err = check_level1_func(ind, fdot, X, Y)
	if err != nil {
		return
	}
	if ind.Nx == 0 {
		return
	}
	Xa := X.ComplexArray()
	Ya := Y.ComplexArray()
	v = zdotu(ind.Nx, Xa[ind.OffsetX:], ind.IncX, Ya[ind.OffsetY:], ind.IncY)
	return
}

// Calculate conjugate X.T * Y for complex matrix.

// See function Dotc.
func DotcComplex(X, Y *matrix.ComplexMatrix, opts ...linalg.Option) (v complex128, err error) {
	v = 0.0
	ind := linalg.GetIndexOpts(opts...)
	err = check_level1_func(ind, fdot, X, Y)
	if err != nil {
		return
	}
	if ind.Nx == 0 {
		return
	}
	Xa := X.ComplexArray()
	Ya := Y.ComplexArray()
	v = zdotc(ind.Nx, Xa[ind.OffsetX:], ind.IncX, Ya[ind.OffsetY:], ind.IncY)
	return
}

// Local Variables:
// tab-width: 4
// End:
