// Copyright (c) Harri Rautila, 2012,2013

// This file is part of github.com/hrautila/linalg/blas package. 
// It is free software, distributed under the terms of GNU Lesser General Public 
// License Version 3, or any later version. See the COPYING tile included in this archive.

package blas

// #cgo LDFLAGS: -L/usr/lib/libblas -lblas
// #include <stdlib.h>
// #include "blas.h"
import "C"
import "unsafe"

//import L "linalg"

// ===========================================================================
// BLAS level 1
// ===========================================================================
// vector - vector 

// Calculate norm2(X). 
func dnrm2(N int, X []float64, incX int) float64 {
    var val C.double
    val = C.dnrm2_((*C.int)(unsafe.Pointer(&N)),
        (*C.double)(unsafe.Pointer(&X[0])),
        (*C.int)(unsafe.Pointer(&incX)))
    return float64(val)
}

// Calculates asum(X). 
func dasum(N int, X []float64, incX int) float64 {
    var val C.double
    val = C.dasum_((*C.int)(unsafe.Pointer(&N)),
        (*C.double)(unsafe.Pointer(&X[0])),
        (*C.int)(unsafe.Pointer(&incX)))
    return float64(val)
}

// Calculates X.T*Y
func ddot(N int, X []float64, incX int, Y []float64, incY int) float64 {
    var val C.double
    val = C.ddot_((*C.int)(unsafe.Pointer(&N)),
        (*C.double)(unsafe.Pointer(&X[0])),
        (*C.int)(unsafe.Pointer(&incX)),
        (*C.double)(unsafe.Pointer(&Y[0])),
        (*C.int)(unsafe.Pointer(&incY)))
    return float64(val)
}

// Find MAX(X) and return index to it.
func idamax(N int, X []float64, incX int) int {
    var idx C.int
    idx = C.idamax_((*C.int)(unsafe.Pointer(&N)),
        (*C.double)(unsafe.Pointer(&X[0])),
        (*C.int)(unsafe.Pointer(&incX)))
    return int(idx)
}

// Swap Y <=> X
func dswap(N int, X []float64, incX int, Y []float64, incY int) {
    C.dswap_((*C.int)(unsafe.Pointer(&N)),
        (*C.double)(unsafe.Pointer(&X[0])),
        (*C.int)(unsafe.Pointer(&incX)),
        (*C.double)(unsafe.Pointer(&Y[0])),
        (*C.int)(unsafe.Pointer(&incY)))
}

// Copy Y <= X
func dcopy(N int, X []float64, incX int, Y []float64, incY int) {
    C.dcopy_((*C.int)(unsafe.Pointer(&N)),
        (*C.double)(unsafe.Pointer(&X[0])),
        (*C.int)(unsafe.Pointer(&incX)),
        (*C.double)(unsafe.Pointer(&Y[0])),
        (*C.int)(unsafe.Pointer(&incY)))
}

// Calculates Y = alpha*X + Y.
func daxpy(N int, alpha float64, X []float64, incX int, Y []float64, incY int) {

    C.daxpy_((*C.int)(unsafe.Pointer(&N)),
        (*C.double)(unsafe.Pointer(&alpha)),
        (*C.double)(unsafe.Pointer(&X[0])),
        (*C.int)(unsafe.Pointer(&incX)),
        (*C.double)(unsafe.Pointer(&Y[0])),
        (*C.int)(unsafe.Pointer(&incY)))
}

// Calculate X = alpha*X.
func dscal(N int, alpha float64, X []float64, incX int) {
    C.dscal_((*C.int)(unsafe.Pointer(&N)),
        (*C.double)(unsafe.Pointer(&alpha)),
        (*C.double)(unsafe.Pointer(&X[0])),
        (*C.int)(unsafe.Pointer(&incX)))
}

/* ------------------------------------------------------------------
 * left out for the time being ....

func drotg(a, b, c, d *float64) {
	C.drotg_((*C.double)(unsafe.Pointer(a)),
		(*C.double)(unsafe.Pointer(b)),
		(*C.double)(unsafe.Pointer(c)),
		(*C.double)(unsafe.Pointer(d)))
}

func drotmg(d1, d2, b1 *float64, b2 float64, P []float64) {
	C.drotmg_((*C.double)(unsafe.Pointer(d1)),
		(*C.double)(unsafe.Pointer(d2)),
		(*C.double)(unsafe.Pointer(b1)),
		(*C.double)(unsafe.Pointer(&b2)),
		(*C.double)(unsafe.Pointer(&P[0])))
}

func drot(N int, X []float64, incX int, Y []float64, incY int, c, s float64) {
	C.drot_((*C.int)(unsafe.Pointer(&N)),
		(*C.double)(unsafe.Pointer(&X[0])),
		(*C.int)(unsafe.Pointer(&incX)),
		(*C.double)(unsafe.Pointer(&Y[0])),
		(*C.int)(unsafe.Pointer(&incY)),
		(*C.double)(unsafe.Pointer(&c)),
		(*C.double)(unsafe.Pointer(&s)))
}

func drotm(N int, X []float64, incX int, Y []float64, incY int, P []float64) {
	C.drotm_((*C.int)(unsafe.Pointer(&N)),
		(*C.double)(unsafe.Pointer(&X[0])),
		(*C.int)(unsafe.Pointer(&incX)),
		(*C.double)(unsafe.Pointer(&Y[0])),
		(*C.int)(unsafe.Pointer(&incY)),
		(*C.double)(unsafe.Pointer(&Y[0])))
}
*/

// ===========================================================================
// BLAS level 2
// ===========================================================================
// Matrix - vector 

// For general matrix A and vector X and Y compute
// Y = alpha * A * X + beta * Y or Y = alpha * A.T * X + beta * Y
func dgemv(transA string, M int, N int, alpha float64,
    A []float64, lda int, X []float64, incX int, beta float64,
    Y []float64, incY int) {

    ctransA := C.CString(transA)
    defer C.free(unsafe.Pointer(ctransA))

    // ?? TODO: protect against index out of bounds panics. 
    C.dgemv_(ctransA,
        (*C.int)(unsafe.Pointer(&M)),
        (*C.int)(unsafe.Pointer(&N)),
        (*C.double)(unsafe.Pointer(&alpha)),
        (*C.double)(unsafe.Pointer(&A[0])),
        (*C.int)(unsafe.Pointer(&lda)),
        (*C.double)(unsafe.Pointer(&X[0])),
        (*C.int)(unsafe.Pointer(&incX)),
        (*C.double)(unsafe.Pointer(&beta)),
        (*C.double)(unsafe.Pointer(&Y[0])),
        (*C.int)(unsafe.Pointer(&incY)))
}

// For general band matrix A  and vector X and Y compute
// Y = alpha * A * X + beta * Y or Y = alpha * A.T * X + beta * Y
func dgbmv(transA string, M int, N int, KL int, KU int,
    alpha float64, A []float64, lda int,
    X []float64, incX int, beta float64,
    Y []float64, incY int) {

    ctransA := C.CString(transA)
    defer C.free(unsafe.Pointer(ctransA))

    // ?? TODO: protect against index out of bounds panics. 
    C.dgbmv_(ctransA,
        (*C.int)(unsafe.Pointer(&M)),
        (*C.int)(unsafe.Pointer(&N)),
        (*C.int)(unsafe.Pointer(&KU)),
        (*C.int)(unsafe.Pointer(&KL)),
        (*C.double)(unsafe.Pointer(&alpha)),
        (*C.double)(unsafe.Pointer(&A[0])),
        (*C.int)(unsafe.Pointer(&lda)),
        (*C.double)(unsafe.Pointer(&X[0])),
        (*C.int)(unsafe.Pointer(&incX)),
        (*C.double)(unsafe.Pointer(&beta)),
        (*C.double)(unsafe.Pointer(&Y[0])),
        (*C.int)(unsafe.Pointer(&incY)))
}

// For triangular matrix A and vector X compute
// X = A * X, X = A.T * X
func dtrmv(uplo, transA, diag string,
    N int, A []float64, lda int, X []float64, incX int) {

    cuplo := C.CString(uplo)
    defer C.free(unsafe.Pointer(cuplo))
    ctransA := C.CString(transA)
    defer C.free(unsafe.Pointer(ctransA))
    cdiag := C.CString(diag)
    defer C.free(unsafe.Pointer(cdiag))

    // ?? TODO: protect against index out of bounds panics. 
    C.dtrmv_(cuplo, ctransA, cdiag,
        (*C.int)(unsafe.Pointer(&N)),
        (*C.double)(unsafe.Pointer(&A[0])),
        (*C.int)(unsafe.Pointer(&lda)),
        (*C.double)(unsafe.Pointer(&X[0])),
        (*C.int)(unsafe.Pointer(&incX)))
}

// For triangular band matrix A and vector X compute
// X = A * X, X = A.T * X
func dtbmv(uplo, transA, diag string,
    N int, K int, A []float64, lda int, X []float64, incX int) {

    cuplo := C.CString(uplo)
    defer C.free(unsafe.Pointer(cuplo))
    ctransA := C.CString(transA)
    defer C.free(unsafe.Pointer(ctransA))
    cdiag := C.CString(diag)
    defer C.free(unsafe.Pointer(cdiag))

    // ?? TODO: protect against index out of bounds panics. 
    C.dtbmv_(cuplo, ctransA, cdiag,
        (*C.int)(unsafe.Pointer(&N)),
        (*C.int)(unsafe.Pointer(&K)),
        (*C.double)(unsafe.Pointer(&A[0])),
        (*C.int)(unsafe.Pointer(&lda)),
        (*C.double)(unsafe.Pointer(&X[0])),
        (*C.int)(unsafe.Pointer(&incX)))
}

// For triangular packed matrix A and vector X compute
// X = A * X, X = A.T * X
func dtpmv(uplo, transA, diag string,
    N int, Ap []float64, X []float64, incX int) {

    cuplo := C.CString(uplo)
    defer C.free(unsafe.Pointer(cuplo))
    ctransA := C.CString(transA)
    defer C.free(unsafe.Pointer(ctransA))
    cdiag := C.CString(diag)
    defer C.free(unsafe.Pointer(cdiag))

    // ?? TODO: protect against index out of bounds panics. 
    C.dtpmv_(cuplo, ctransA, cdiag,
        (*C.int)(unsafe.Pointer(&N)),
        (*C.double)(unsafe.Pointer(&Ap[0])),
        (*C.double)(unsafe.Pointer(&X[0])),
        (*C.int)(unsafe.Pointer(&incX)))
}

// For triangular matrix A and vector X solve
// X = inv(A) * X or X = inv(A.T) * X
func dtrsv(uplo, transA, diag string,
    N int, A []float64, lda int, X []float64, incX int) {

    cuplo := C.CString(uplo)
    defer C.free(unsafe.Pointer(cuplo))
    ctransA := C.CString(transA)
    defer C.free(unsafe.Pointer(ctransA))
    cdiag := C.CString(diag)
    defer C.free(unsafe.Pointer(cdiag))

    // ?? TODO: protect against index out of bounds panics. 
    C.dtrsv_(cuplo, ctransA, cdiag,
        (*C.int)(unsafe.Pointer(&N)),
        (*C.double)(unsafe.Pointer(&A[0])),
        (*C.int)(unsafe.Pointer(&lda)),
        (*C.double)(unsafe.Pointer(&X[0])),
        (*C.int)(unsafe.Pointer(&incX)))
}

// For triangular band matrix A and vector X solve
// X = inv(A) * X or X = inv(A.T) * X
func dtbsv(uplo, transA, diag string,
    N int, K int, A []float64, lda int, X []float64, incX int) {

    cuplo := C.CString(uplo)
    defer C.free(unsafe.Pointer(cuplo))
    ctransA := C.CString(transA)
    defer C.free(unsafe.Pointer(ctransA))
    cdiag := C.CString(diag)
    defer C.free(unsafe.Pointer(cdiag))

    // ?? TODO: protect against index out of bounds panics. 
    C.dtbsv_(cuplo, ctransA, cdiag,
        (*C.int)(unsafe.Pointer(&N)),
        (*C.int)(unsafe.Pointer(&K)),
        (*C.double)(unsafe.Pointer(&A[0])),
        (*C.int)(unsafe.Pointer(&lda)),
        (*C.double)(unsafe.Pointer(&X[0])),
        (*C.int)(unsafe.Pointer(&incX)))
}

// For triangular matrix A and vector X solve
// X = inv(A) * X or X = inv(A.T) * X
func dtpsv(uplo, transA, diag string,
    N int, Ap []float64, X []float64, incX int) {

    cuplo := C.CString(uplo)
    defer C.free(unsafe.Pointer(cuplo))
    ctransA := C.CString(transA)
    defer C.free(unsafe.Pointer(ctransA))
    cdiag := C.CString(diag)
    defer C.free(unsafe.Pointer(cdiag))

    // ?? TODO: protect against index out of bounds panics. 
    C.dtpsv_(cuplo, ctransA, cdiag,
        (*C.int)(unsafe.Pointer(&N)),
        (*C.double)(unsafe.Pointer(&Ap[0])),
        (*C.double)(unsafe.Pointer(&X[0])),
        (*C.int)(unsafe.Pointer(&incX)))
}

// For symmetric matrix  and vector X solve
// Y = alpha * A * X + beta * Y
func dsymv(uplo string, N int, alpha float64,
    A []float64, lda int, X []float64, incX int, beta float64,
    Y []float64, incY int) {

    cuplo := C.CString(uplo)
    defer C.free(unsafe.Pointer(cuplo))

    // ?? TODO: protect against index out of bounds panics. 
    C.dsymv_(cuplo, (*C.int)(unsafe.Pointer(&N)),
        (*C.double)(unsafe.Pointer(&alpha)),
        (*C.double)(unsafe.Pointer(&A[0])),
        (*C.int)(unsafe.Pointer(&lda)),
        (*C.double)(unsafe.Pointer(&X[0])),
        (*C.int)(unsafe.Pointer(&incX)),
        (*C.double)(unsafe.Pointer(&beta)),
        (*C.double)(unsafe.Pointer(&Y[0])),
        (*C.int)(unsafe.Pointer(&incY)))

}

// For symmetric band matrix  and vector X solve
// Y = alpha * A * X + beta * Y
func dsbmv(uplo string, N int, K int, alpha float64,
    A []float64, lda int, X []float64, incX int, beta float64,
    Y []float64, incY int) {

    cuplo := C.CString(uplo)
    defer C.free(unsafe.Pointer(cuplo))

    // ?? TODO: protect against index out of bounds panics. 
    C.dsbmv_(cuplo,
        (*C.int)(unsafe.Pointer(&N)),
        (*C.int)(unsafe.Pointer(&K)),
        (*C.double)(unsafe.Pointer(&alpha)),
        (*C.double)(unsafe.Pointer(&A[0])),
        (*C.int)(unsafe.Pointer(&lda)),
        (*C.double)(unsafe.Pointer(&X[0])),
        (*C.int)(unsafe.Pointer(&incX)),
        (*C.double)(unsafe.Pointer(&beta)),
        (*C.double)(unsafe.Pointer(&Y[0])),
        (*C.int)(unsafe.Pointer(&incY)))

}

// For symmetric packed matrix  and vector X solve
// Y = alpha * A * X + beta * Y

func dspmv(uplo string, N int, alpha float64,
    Ap []float64, X []float64, incX int, beta float64,
    Y []float64, incY int) {

    cuplo := C.CString(uplo)
    defer C.free(unsafe.Pointer(cuplo))

    // ?? TODO: protect against index out of bounds panics. 
    C.dspmv_(cuplo,
        (*C.int)(unsafe.Pointer(&N)),
        (*C.double)(unsafe.Pointer(&alpha)),
        (*C.double)(unsafe.Pointer(&Ap[0])),
        (*C.double)(unsafe.Pointer(&X[0])),
        (*C.int)(unsafe.Pointer(&incX)),
        (*C.double)(unsafe.Pointer(&beta)),
        (*C.double)(unsafe.Pointer(&Y[0])),
        (*C.int)(unsafe.Pointer(&incY)))

}

func dger(M int, N int, alpha float64,
    X []float64, incX int, Y []float64, incY int,
    A []float64, lda int) {

    // ?? TODO: protect against index out of bounds panics. 
    C.dger_((*C.int)(unsafe.Pointer(&M)),
        (*C.int)(unsafe.Pointer(&N)),
        (*C.double)(unsafe.Pointer(&alpha)),
        (*C.double)(unsafe.Pointer(&X[0])),
        (*C.int)(unsafe.Pointer(&incX)),
        (*C.double)(unsafe.Pointer(&Y[0])),
        (*C.int)(unsafe.Pointer(&incY)),
        (*C.double)(unsafe.Pointer(&A[0])),
        (*C.int)(unsafe.Pointer(&lda)))

}

func dsyr(uplo string, N int, alpha float64,
    X []float64, incX int, A []float64, lda int) {

    cuplo := C.CString(uplo)
    defer C.free(unsafe.Pointer(cuplo))

    // ?? TODO: protect against index out of bounds panics. 
    C.dsyr_(cuplo,
        (*C.int)(unsafe.Pointer(&N)),
        (*C.double)(unsafe.Pointer(&alpha)),
        (*C.double)(unsafe.Pointer(&X[0])),
        (*C.int)(unsafe.Pointer(&incX)),
        (*C.double)(unsafe.Pointer(&A[0])),
        (*C.int)(unsafe.Pointer(&lda)))

}

func dspr(uplo string, N int, alpha float64,
    X []float64, incX int, Ap []float64) {

    cuplo := C.CString(uplo)
    defer C.free(unsafe.Pointer(cuplo))

    // ?? TODO: protect against index out of bounds panics. 
    C.dspr_(cuplo,
        (*C.int)(unsafe.Pointer(&N)),
        (*C.double)(unsafe.Pointer(&alpha)),
        (*C.double)(unsafe.Pointer(&X[0])),
        (*C.int)(unsafe.Pointer(&incX)),
        (*C.double)(unsafe.Pointer(&Ap[0])))

}

func dsyr2(uplo string, N int, alpha float64,
    X []float64, incX int, Y []float64, incY int, A []float64, lda int) {

    cuplo := C.CString(uplo)
    defer C.free(unsafe.Pointer(cuplo))

    // ?? TODO: protect against index out of bounds panics. 

    C.dsyr2_(cuplo,
        (*C.int)(unsafe.Pointer(&N)),
        (*C.double)(unsafe.Pointer(&alpha)),
        (*C.double)(unsafe.Pointer(&X[0])),
        (*C.int)(unsafe.Pointer(&incX)),
        (*C.double)(unsafe.Pointer(&Y[0])),
        (*C.int)(unsafe.Pointer(&incY)),
        (*C.double)(unsafe.Pointer(&A[0])),
        (*C.int)(unsafe.Pointer(&lda)))

}

func dspr2(uplo string, N int, alpha float64,
    X []float64, incX int, Y []float64, incY int, Ap []float64) {

    cuplo := C.CString(uplo)
    defer C.free(unsafe.Pointer(cuplo))

    // ?? TODO: protect against index out of bounds panics. 
    C.dspr2_(cuplo,
        (*C.int)(unsafe.Pointer(&N)),
        (*C.double)(unsafe.Pointer(&alpha)),
        (*C.double)(unsafe.Pointer(&X[0])),
        (*C.int)(unsafe.Pointer(&incX)),
        (*C.double)(unsafe.Pointer(&Y[0])),
        (*C.int)(unsafe.Pointer(&incY)),
        (*C.double)(unsafe.Pointer(&Ap[0])))

}

// ===========================================================================
// BLAS level 3
// ===========================================================================
// Matrix - Matrix

func dgemm(transA, transB string, M int, N int, K int,
    alpha float64, A []float64, lda int, B []float64, ldb int, beta float64,
    C []float64, ldc int) {

    ctransA := C.CString(transA)
    defer C.free(unsafe.Pointer(ctransA))
    ctransB := C.CString(transB)
    defer C.free(unsafe.Pointer(ctransB))

    // protect against index out of bounds panics
    var aptr, bptr, cptr *float64 = nil, nil, nil
    if len(A) > 0 {
        aptr = &A[0]
    }
    if len(B) > 0 {
        bptr = &B[0]
    }
    if len(C) > 0 {
        cptr = &C[0]
    }
    C.dgemm_(ctransA, ctransB,
        (*C.int)(unsafe.Pointer(&M)),
        (*C.int)(unsafe.Pointer(&N)),
        (*C.int)(unsafe.Pointer(&K)),
        (*C.double)(unsafe.Pointer(&alpha)),
        (*C.double)(unsafe.Pointer(aptr)),
        (*C.int)(unsafe.Pointer(&lda)),
        (*C.double)(unsafe.Pointer(bptr)),
        (*C.int)(unsafe.Pointer(&ldb)),
        (*C.double)(unsafe.Pointer(&beta)),
        (*C.double)(unsafe.Pointer(cptr)),
        (*C.int)(unsafe.Pointer(&ldc)))
}

func dsymm(side, uplo string, M int, N int,
    alpha float64, A []float64, lda int, B []float64, ldb int, beta float64,
    C []float64, ldc int) {

    cuplo := C.CString(uplo)
    defer C.free(unsafe.Pointer(cuplo))
    cside := C.CString(side)
    defer C.free(unsafe.Pointer(cside))

    // protect against index out of bounds panics
    var aptr, bptr, cptr *float64 = nil, nil, nil
    if len(A) > 0 {
        aptr = &A[0]
    }
    if len(B) > 0 {
        bptr = &B[0]
    }
    if len(C) > 0 {
        cptr = &C[0]
    }
    C.dsymm_(cside, cuplo,
        (*C.int)(unsafe.Pointer(&M)),
        (*C.int)(unsafe.Pointer(&N)),
        (*C.double)(unsafe.Pointer(&alpha)),
        (*C.double)(unsafe.Pointer(aptr)),
        (*C.int)(unsafe.Pointer(&lda)),
        (*C.double)(unsafe.Pointer(bptr)),
        (*C.int)(unsafe.Pointer(&ldb)),
        (*C.double)(unsafe.Pointer(&beta)),
        (*C.double)(unsafe.Pointer(cptr)),
        (*C.int)(unsafe.Pointer(&ldc)))
}

func dsyrk(uplo, trans string, N int, K int,
    alpha float64, A []float64, lda int, beta float64,
    C []float64, ldc int) {

    cuplo := C.CString(uplo)
    defer C.free(unsafe.Pointer(cuplo))
    ctrans := C.CString(trans)
    defer C.free(unsafe.Pointer(ctrans))

    // protect against index out of bounds panics
    var aptr, cptr *float64 = nil, nil
    if len(A) > 0 {
        aptr = &A[0]
    }
    if len(C) > 0 {
        cptr = &C[0]
    }
    C.dsyrk_(cuplo, ctrans,
        (*C.int)(unsafe.Pointer(&N)),
        (*C.int)(unsafe.Pointer(&K)),
        (*C.double)(unsafe.Pointer(&alpha)),
        (*C.double)(unsafe.Pointer(aptr)),
        (*C.int)(unsafe.Pointer(&lda)),
        (*C.double)(unsafe.Pointer(&beta)),
        (*C.double)(unsafe.Pointer(cptr)),
        (*C.int)(unsafe.Pointer(&ldc)))
}

func dsyr2k(uplo, trans string, N int, K int,
    alpha float64, A []float64, lda int, B []float64, ldb int, beta float64,
    C []float64, ldc int) {

    cuplo := C.CString(uplo)
    defer C.free(unsafe.Pointer(cuplo))
    ctrans := C.CString(trans)
    defer C.free(unsafe.Pointer(ctrans))

    // protect against index out of bounds panics
    var aptr, bptr, cptr *float64 = nil, nil, nil
    if len(A) > 0 {
        aptr = &A[0]
    }
    if len(B) > 0 {
        bptr = &B[0]
    }
    if len(C) > 0 {
        cptr = &C[0]
    }
    C.dsyr2k_(cuplo, ctrans,
        (*C.int)(unsafe.Pointer(&N)),
        (*C.int)(unsafe.Pointer(&K)),
        (*C.double)(unsafe.Pointer(&alpha)),
        (*C.double)(unsafe.Pointer(aptr)),
        (*C.int)(unsafe.Pointer(&lda)),
        (*C.double)(unsafe.Pointer(bptr)),
        (*C.int)(unsafe.Pointer(&ldb)),
        (*C.double)(unsafe.Pointer(&beta)),
        (*C.double)(unsafe.Pointer(cptr)),
        (*C.int)(unsafe.Pointer(&ldc)))
}

func dtrmm(side, uplo, transA, diag string,
    M int, N int, alpha float64, A []float64, lda int, B []float64, ldb int) {

    cuplo := C.CString(uplo)
    defer C.free(unsafe.Pointer(cuplo))
    ctrans := C.CString(transA)
    defer C.free(unsafe.Pointer(ctrans))
    cside := C.CString(side)
    defer C.free(unsafe.Pointer(cside))
    cdiag := C.CString(diag)
    defer C.free(unsafe.Pointer(cdiag))

    // protect against index out of bounds panics
    var aptr, bptr *float64 = nil, nil
    if len(A) > 0 {
        aptr = &A[0]
    }
    if len(B) > 0 {
        bptr = &B[0]
    }
    C.dtrmm_(cside, cuplo, ctrans, cdiag,
        (*C.int)(unsafe.Pointer(&M)),
        (*C.int)(unsafe.Pointer(&N)),
        (*C.double)(unsafe.Pointer(&alpha)),
        (*C.double)(unsafe.Pointer(aptr)),
        (*C.int)(unsafe.Pointer(&lda)),
        (*C.double)(unsafe.Pointer(bptr)),
        (*C.int)(unsafe.Pointer(&ldb)))
}

func dtrsm(side, uplo, transA, diag string,
    M int, N int, alpha float64, A []float64, lda int, B []float64, ldb int) {

    cuplo := C.CString(uplo)
    defer C.free(unsafe.Pointer(cuplo))
    ctrans := C.CString(transA)
    defer C.free(unsafe.Pointer(ctrans))
    cside := C.CString(side)
    defer C.free(unsafe.Pointer(cside))
    cdiag := C.CString(diag)
    defer C.free(unsafe.Pointer(cdiag))

    // protect against index out of bounds panics
    var aptr, bptr *float64 = nil, nil
    if len(A) > 0 {
        aptr = &A[0]
    }
    if len(B) > 0 {
        bptr = &B[0]
    }
    C.dtrsm_(cside, cuplo, ctrans, cdiag,
        (*C.int)(unsafe.Pointer(&M)),
        (*C.int)(unsafe.Pointer(&N)),
        (*C.double)(unsafe.Pointer(&alpha)),
        (*C.double)(unsafe.Pointer(aptr)),
        (*C.int)(unsafe.Pointer(&lda)),
        (*C.double)(unsafe.Pointer(bptr)),
        (*C.int)(unsafe.Pointer(&ldb)))
}

// Local Variables:
// tab-width: 4
// End:
