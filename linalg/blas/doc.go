// Copyright (c) Harri Rautila, 2012,2013

// This file is part of go.opt/linalg package. It is free software, distributed
// under the terms of GNU Lesser General Public License Version 3, or any later
// version. See the COPYING tile included in this archive.

// Interface to the double-precision real and complex BLAS library.
//
// This package is implementation of CVXOPT Python blas interface in GO.
//
// Double and complex matrices and vectors are stored in column major 
// matrices using the conventional BLAS storage schemes, with the
// matrix buffers interpreted as one-dimensional arrays.
// 
// For each matrix argument X, an additional integer argument
// offsetX specifies the start of the array, i.e., the pointer
// of X[offsetX:] is passed to the BLAS function.  The other 
// arguments (dimensions and options) have the same meaning as in
// the BLAS definition.  Default values of the dimension arguments
// are derived from the matrix sizes.
package blas
