// Copyright (c) Harri Rautila, 2012,2013

// This file is part of github.com/hrautila/linalg/blas package. 
// It is free software, distributed under the terms of GNU Lesser General Public 
// License Version 3, or any later version. See the COPYING tile included in this archive.

package blas

// #cgo LDFLAGS: -L/usr/lib/libblas -lblas -L/usr/lib/atlas -latlas
// #include <stdlib.h>
// #include "blas.h"
import "C"
import "unsafe"

// ===========================================================================
// BLAS level 1
// ===========================================================================
// vector - vector 

// Calculate norm2(X). 
func dznrm2(N int, X []complex128, incX int) float64 {
    var val C.double
    val = C.dznrm2_((*C.int)(unsafe.Pointer(&N)),
        unsafe.Pointer(&X[0]),
        (*C.int)(unsafe.Pointer(&incX)))
    return float64(val)
}

// Calculates asum(X). 
func dzasum(N int, X []complex128, incX int) float64 {
    var val C.double
    val = C.dzasum_((*C.int)(unsafe.Pointer(&N)),
        unsafe.Pointer(&X[0]),
        (*C.int)(unsafe.Pointer(&incX)))
    return float64(val)
}

// Calculates conjucate X.T*Y
func zdotc(N int, X []complex128, incX int, Y []complex128, incY int) complex128 {
    var result complex128
    C.zdotcsub_((*C.int)(unsafe.Pointer(&N)),
        unsafe.Pointer(&X[0]),
        (*C.int)(unsafe.Pointer(&incX)),
        unsafe.Pointer(&Y[0]),
        (*C.int)(unsafe.Pointer(&incY)),
        unsafe.Pointer(&result))
    return result
}

// Calculates X.T*Y
func zdotu(N int, X []complex128, incX int, Y []complex128, incY int) complex128 {
    var result complex128
    C.zdotusub_((*C.int)(unsafe.Pointer(&N)),
        unsafe.Pointer(&X[0]),
        (*C.int)(unsafe.Pointer(&incX)),
        unsafe.Pointer(&Y[0]),
        (*C.int)(unsafe.Pointer(&incY)),
        unsafe.Pointer(&result))
    return result
}

// Find MAX(X) and return index to it.
func izamax(N int, X []complex128, incX int) int {
    var idx C.int
    idx = C.izamax_((*C.int)(unsafe.Pointer(&N)),
        unsafe.Pointer(&X[0]),
        (*C.int)(unsafe.Pointer(&incX)))
    return int(idx)
}

// Swap Y <=> X
func zswap(N int, X []complex128, incX int, Y []complex128, incY int) {
    C.zswap_((*C.int)(unsafe.Pointer(&N)),
        unsafe.Pointer(&X[0]),
        (*C.int)(unsafe.Pointer(&incX)),
        unsafe.Pointer(&Y[0]),
        (*C.int)(unsafe.Pointer(&incY)))
}

// Copy Y <= X
func zcopy(N int, X []complex128, incX int, Y []complex128, incY int) {
    C.zcopy_((*C.int)(unsafe.Pointer(&N)),
        unsafe.Pointer(&X[0]),
        (*C.int)(unsafe.Pointer(&incX)),
        unsafe.Pointer(&Y[0]),
        (*C.int)(unsafe.Pointer(&incY)))
}

// Calculates Y = alpha*X + Y.
func zaxpy(N int, alpha complex128, X []complex128, incX int, Y []complex128, incY int) {

    C.zaxpy_((*C.int)(unsafe.Pointer(&N)),
        unsafe.Pointer(&alpha),
        unsafe.Pointer(&X[0]),
        (*C.int)(unsafe.Pointer(&incX)),
        unsafe.Pointer(&Y[0]),
        (*C.int)(unsafe.Pointer(&incY)))
}

// Calculate X = alpha*X.
func zscal(N int, alpha complex128, X []complex128, incX int) {
    C.zscal_((*C.int)(unsafe.Pointer(&N)),
        (unsafe.Pointer(&alpha)),
        (unsafe.Pointer(&X[0])),
        (*C.int)(unsafe.Pointer(&incX)))
}

/* ------------------------------------------------------------------
 * left out for the time being ....

func zrotg(a, b, c, d *complex128) {
	C.zrotg_((unsafe.Pointer(a)),
		(unsafe.Pointer(b)),
		(unsafe.Pointer(c)),
		(unsafe.Pointer(d)))
}

func zrotmg(d1, d2, b1 *complex128, b2 complex128, P []complex128) {
	C.zrotmg_((unsafe.Pointer(d1)),
		(unsafe.Pointer(d2)),
		(unsafe.Pointer(b1)),
		(unsafe.Pointer(&b2)),
		(unsafe.Pointer(&P[0])))
}

func zrot(N int, X []complex128, incX int, Y []complex128, incY int, c, s complex128) {
	C.zrot_((*C.int)(unsafe.Pointer(&N)),
		(unsafe.Pointer(&X[0])),
		(*C.int)(unsafe.Pointer(&incX)),
		(unsafe.Pointer(&Y[0])),
		(*C.int)(unsafe.Pointer(&incY)),
		(unsafe.Pointer(&c)),
		(unsafe.Pointer(&s)))
}

func zrotm(N int, X []complex128, incX int, Y []complex128, incY int, P []complex128) {
	C.zrotm_((*C.int)(unsafe.Pointer(&N)),
		(unsafe.Pointer(&X[0])),
		(*C.int)(unsafe.Pointer(&incX)),
		(unsafe.Pointer(&Y[0])),
		(*C.int)(unsafe.Pointer(&incY)),
		(unsafe.Pointer(&Y[0])))
}
*/

// ===========================================================================
// BLAS level 2
// ===========================================================================
// Matrix - vector 

// For general matrix A and vector X and Y compute
// Y = alpha * A * X + beta * Y or Y = alpha * A.T * X + beta * Y
func zgemv(transA string, M int, N int, alpha complex128,
    A []complex128, lda int, X []complex128, incX int, beta complex128,
    Y []complex128, incY int) {

    ctransA := C.CString(transA)
    defer C.free(unsafe.Pointer(ctransA))

    C.zgemv_(ctransA,
        (*C.int)(unsafe.Pointer(&M)),
        (*C.int)(unsafe.Pointer(&N)),
        (unsafe.Pointer(&alpha)),
        (unsafe.Pointer(&A[0])),
        (*C.int)(unsafe.Pointer(&lda)),
        (unsafe.Pointer(&X[0])),
        (*C.int)(unsafe.Pointer(&incX)),
        (unsafe.Pointer(&beta)),
        (unsafe.Pointer(&Y[0])),
        (*C.int)(unsafe.Pointer(&incY)))
}

// For general band matrix A  and vector X and Y compute
// Y = alpha * A * X + beta * Y or Y = alpha * A.T * X + beta * Y
func zgbmv(transA string, M int, N int, KL int, KU int,
    alpha complex128, A []complex128, lda int,
    X []complex128, incX int, beta complex128,
    Y []complex128, incY int) {

    ctransA := C.CString(transA)
    defer C.free(unsafe.Pointer(ctransA))

    C.zgbmv_(ctransA,
        (*C.int)(unsafe.Pointer(&M)),
        (*C.int)(unsafe.Pointer(&N)),
        (*C.int)(unsafe.Pointer(&KU)),
        (*C.int)(unsafe.Pointer(&KL)),
        (unsafe.Pointer(&alpha)),
        (unsafe.Pointer(&A[0])),
        (*C.int)(unsafe.Pointer(&lda)),
        (unsafe.Pointer(&X[0])),
        (*C.int)(unsafe.Pointer(&incX)),
        (unsafe.Pointer(&beta)),
        (unsafe.Pointer(&Y[0])),
        (*C.int)(unsafe.Pointer(&incY)))
}

// For triangular matrix A and vector X compute
// X = A * X, X = A.T * X
func ztrmv(uplo, transA, diag string,
    N int, A []complex128, lda int, X []complex128, incX int) {

    cuplo := C.CString(uplo)
    defer C.free(unsafe.Pointer(cuplo))
    ctransA := C.CString(transA)
    defer C.free(unsafe.Pointer(ctransA))
    cdiag := C.CString(diag)
    defer C.free(unsafe.Pointer(cdiag))

    C.ztrmv_(cuplo, ctransA, cdiag,
        (*C.int)(unsafe.Pointer(&N)),
        (unsafe.Pointer(&A[0])),
        (*C.int)(unsafe.Pointer(&lda)),
        (unsafe.Pointer(&X[0])),
        (*C.int)(unsafe.Pointer(&incX)))
}

// For triangular band matrix A and vector X compute
// X = A * X, X = A.T * X
func ztbmv(uplo, transA, diag string,
    N int, K int, A []complex128, lda int, X []complex128, incX int) {

    cuplo := C.CString(uplo)
    defer C.free(unsafe.Pointer(cuplo))
    ctransA := C.CString(transA)
    defer C.free(unsafe.Pointer(ctransA))
    cdiag := C.CString(diag)
    defer C.free(unsafe.Pointer(cdiag))

    C.ztbmv_(cuplo, ctransA, cdiag,
        (*C.int)(unsafe.Pointer(&N)),
        (*C.int)(unsafe.Pointer(&K)),
        (unsafe.Pointer(&A[0])),
        (*C.int)(unsafe.Pointer(&lda)),
        (unsafe.Pointer(&X[0])),
        (*C.int)(unsafe.Pointer(&incX)))
}

// For triangular matrix A and vector X solve
// X = inv(A) * X or X = inv(A.T) * X
func ztrsv(uplo, transA, diag string,
    N int, A []complex128, lda int, X []complex128, incX int) {

    cuplo := C.CString(uplo)
    defer C.free(unsafe.Pointer(cuplo))
    ctransA := C.CString(transA)
    defer C.free(unsafe.Pointer(ctransA))
    cdiag := C.CString(diag)
    defer C.free(unsafe.Pointer(cdiag))

    C.ztrsv_(cuplo, ctransA, cdiag,
        (*C.int)(unsafe.Pointer(&N)),
        (unsafe.Pointer(&A[0])),
        (*C.int)(unsafe.Pointer(&lda)),
        (unsafe.Pointer(&X[0])),
        (*C.int)(unsafe.Pointer(&incX)))
}

// For triangular band matrix A and vector X solve
// X = inv(A) * X or X = inv(A.T) * X
func ztbsv(uplo, transA, diag string,
    N int, K int, A []complex128, lda int, X []complex128, incX int) {

    cuplo := C.CString(uplo)
    defer C.free(unsafe.Pointer(cuplo))
    ctransA := C.CString(transA)
    defer C.free(unsafe.Pointer(ctransA))
    cdiag := C.CString(diag)
    defer C.free(unsafe.Pointer(cdiag))

    C.ztbsv_(cuplo, ctransA, cdiag,
        (*C.int)(unsafe.Pointer(&N)),
        (*C.int)(unsafe.Pointer(&K)),
        (unsafe.Pointer(&A[0])),
        (*C.int)(unsafe.Pointer(&lda)),
        (unsafe.Pointer(&X[0])),
        (*C.int)(unsafe.Pointer(&incX)))
}

// For symmetric matrix  and vector X solve
// Y = alpha * A * X + beta * Y
func zhemv(uplo string, N int, alpha complex128,
    A []complex128, lda int, X []complex128, incX int, beta complex128,
    Y []complex128, incY int) {

    cuplo := C.CString(uplo)
    defer C.free(unsafe.Pointer(cuplo))

    C.zhemv_(cuplo, (*C.int)(unsafe.Pointer(&N)),
        (unsafe.Pointer(&alpha)),
        (unsafe.Pointer(&A[0])),
        (*C.int)(unsafe.Pointer(&lda)),
        (unsafe.Pointer(&X[0])),
        (*C.int)(unsafe.Pointer(&incX)),
        (unsafe.Pointer(&beta)),
        (unsafe.Pointer(&Y[0])),
        (*C.int)(unsafe.Pointer(&incY)))

}

// For symmetric band matrix  and vector X solve
// Y = alpha * A * X + beta * Y
func zhbmv(uplo string, N int, K int, alpha complex128,
    A []complex128, lda int, X []complex128, incX int, beta complex128,
    Y []complex128, incY int) {

    cuplo := C.CString(uplo)
    defer C.free(unsafe.Pointer(cuplo))

    C.zhbmv_(cuplo,
        (*C.int)(unsafe.Pointer(&N)),
        (*C.int)(unsafe.Pointer(&K)),
        (unsafe.Pointer(&alpha)),
        (unsafe.Pointer(&A[0])),
        (*C.int)(unsafe.Pointer(&lda)),
        (unsafe.Pointer(&X[0])),
        (*C.int)(unsafe.Pointer(&incX)),
        (unsafe.Pointer(&beta)),
        (unsafe.Pointer(&Y[0])),
        (*C.int)(unsafe.Pointer(&incY)))

}

func zgerc(M int, N int, alpha complex128,
    X []complex128, incX int, Y []complex128, incY int,
    A []complex128, lda int) {

    C.zgerc_((*C.int)(unsafe.Pointer(&M)),
        (*C.int)(unsafe.Pointer(&N)),
        (unsafe.Pointer(&alpha)),
        (unsafe.Pointer(&X[0])),
        (*C.int)(unsafe.Pointer(&incX)),
        (unsafe.Pointer(&Y[0])),
        (*C.int)(unsafe.Pointer(&incY)),
        (unsafe.Pointer(&A[0])),
        (*C.int)(unsafe.Pointer(&lda)))

}

func zgeru(M int, N int, alpha complex128,
    X []complex128, incX int, Y []complex128, incY int,
    A []complex128, lda int) {

    C.zgeru_((*C.int)(unsafe.Pointer(&M)),
        (*C.int)(unsafe.Pointer(&N)),
        (unsafe.Pointer(&alpha)),
        (unsafe.Pointer(&X[0])),
        (*C.int)(unsafe.Pointer(&incX)),
        (unsafe.Pointer(&Y[0])),
        (*C.int)(unsafe.Pointer(&incY)),
        (unsafe.Pointer(&A[0])),
        (*C.int)(unsafe.Pointer(&lda)))

}

func zher(uplo string, N int, alpha float64,
    X []complex128, incX int, A []complex128, lda int) {

    cuplo := C.CString(uplo)
    defer C.free(unsafe.Pointer(cuplo))

    C.zher_(cuplo,
        (*C.int)(unsafe.Pointer(&N)),
        (*C.double)(unsafe.Pointer(&alpha)),
        (unsafe.Pointer(&X[0])),
        (*C.int)(unsafe.Pointer(&incX)),
        (unsafe.Pointer(&A[0])),
        (*C.int)(unsafe.Pointer(&lda)))

}

func zher2(uplo string, N int, alpha complex128,
    X []complex128, incX int, Y []complex128, incY int, A []complex128, lda int) {

    cuplo := C.CString(uplo)
    defer C.free(unsafe.Pointer(cuplo))

    C.zher2_(cuplo,
        (*C.int)(unsafe.Pointer(&N)),
        (unsafe.Pointer(&alpha)),
        (unsafe.Pointer(&X[0])),
        (*C.int)(unsafe.Pointer(&incX)),
        (unsafe.Pointer(&Y[0])),
        (*C.int)(unsafe.Pointer(&incY)),
        (unsafe.Pointer(&A[0])),
        (*C.int)(unsafe.Pointer(&lda)))

}

/**
// For triangular packed matrix A and vector X compute
// X = A * X, X = A.T * X
func ztpmv(uplo, transA, diag string,
		N int, Ap []complex128, X []complex128, incX int) {

	cuplo := C.CString(uplo)
	defer C.free(unsafe.Pointer(cuplo))
	ctransA := C.CString(transA)
	defer C.free(unsafe.Pointer(ctransA))
	cdiag := C.CString(diag)
	defer C.free(unsafe.Pointer(cdiag))

	C.ztpmv_(cuplo, ctransA, cdiag,
		(*C.int)(unsafe.Pointer(&N)),
		(unsafe.Pointer(&Ap[0])), 
		(unsafe.Pointer(&X[0])),
		(*C.int)(unsafe.Pointer(&incX)))
}


// For triangular matrix A and vector X solve
// X = inv(A) * X or X = inv(A.T) * X
func ztpsv(uplo, transA, diag string,
		N int, Ap []complex128, X []complex128, incX int) {

	cuplo := C.CString(uplo)
	defer C.free(unsafe.Pointer(cuplo))
	ctransA := C.CString(transA)
	defer C.free(unsafe.Pointer(ctransA))
	cdiag := C.CString(diag)
	defer C.free(unsafe.Pointer(cdiag))

	C.ztpsv_(cuplo, ctransA, cdiag,
		(*C.int)(unsafe.Pointer(&N)),
		(unsafe.Pointer(&Ap[0])), 
		(unsafe.Pointer(&X[0])),
		(*C.int)(unsafe.Pointer(&incX)))
}

// For symmetric packed matrix  and vector X solve
// Y = alpha * A * X + beta * Y
func zspmv(uplo string, N int, alpha complex128,
	Ap []complex128, X []complex128, incX int, beta complex128,
	Y []complex128, incY int) {

	cuplo := C.CString(uplo)
	defer C.free(unsafe.Pointer(cuplo))

	C.zspmv_(cuplo,
		(*C.int)(unsafe.Pointer(&N)),
		(unsafe.Pointer(&alpha)),
		(unsafe.Pointer(&Ap[0])), 
		(unsafe.Pointer(&X[0])),
		(*C.int)(unsafe.Pointer(&incX)),
		(unsafe.Pointer(&beta)),
		(unsafe.Pointer(&Y[0])),
		(*C.int)(unsafe.Pointer(&incY)))

}

func zspr(uplo string, N int, alpha complex128,
		X []complex128, incX int, Ap []complex128) {

	cuplo := C.CString(uplo)
	defer C.free(unsafe.Pointer(cuplo))

	C.zspr_(cuplo,
		(*C.int)(unsafe.Pointer(&N)),
		(unsafe.Pointer(&alpha)),
		(unsafe.Pointer(&X[0])),
		(*C.int)(unsafe.Pointer(&incX)),
		(unsafe.Pointer(&Ap[0])))

}

func zspr2(uplo string, N int, alpha complex128,
	X []complex128, incX int, Y []complex128, incY int, Ap []complex128) {

	cuplo := C.CString(uplo)
	defer C.free(unsafe.Pointer(cuplo))

	C.zspr2_(cuplo,
		(*C.int)(unsafe.Pointer(&N)),
		(unsafe.Pointer(&alpha)),
		(unsafe.Pointer(&X[0])),
		(*C.int)(unsafe.Pointer(&incX)),
		(unsafe.Pointer(&Y[0])),
		(*C.int)(unsafe.Pointer(&incY)),
		(unsafe.Pointer(&Ap[0])))

}
**/

// ===========================================================================
// BLAS level 3
// ===========================================================================
// Matrix - Matrix

func zgemm(transA, transB string, M int, N int, K int,
    alpha complex128, A []complex128, lda int, B []complex128, ldb int, beta complex128,
    C []complex128, ldc int) {

    ctransA := C.CString(transA)
    defer C.free(unsafe.Pointer(ctransA))
    ctransB := C.CString(transB)
    defer C.free(unsafe.Pointer(ctransB))

    C.zgemm_(ctransA, ctransB,
        (*C.int)(unsafe.Pointer(&M)),
        (*C.int)(unsafe.Pointer(&N)),
        (*C.int)(unsafe.Pointer(&K)),
        (unsafe.Pointer(&alpha)),
        (unsafe.Pointer(&A[0])),
        (*C.int)(unsafe.Pointer(&lda)),
        (unsafe.Pointer(&B[0])),
        (*C.int)(unsafe.Pointer(&ldb)),
        (unsafe.Pointer(&beta)),
        (unsafe.Pointer(&C[0])),
        (*C.int)(unsafe.Pointer(&ldc)))
}

func zsymm(side, uplo string, M int, N int,
    alpha complex128, A []complex128, lda int, B []complex128, ldb int, beta complex128,
    C []complex128, ldc int) {

    cuplo := C.CString(uplo)
    defer C.free(unsafe.Pointer(cuplo))
    cside := C.CString(side)
    defer C.free(unsafe.Pointer(cside))

    C.zsymm_(cside, cuplo,
        (*C.int)(unsafe.Pointer(&M)),
        (*C.int)(unsafe.Pointer(&N)),
        (unsafe.Pointer(&alpha)),
        (unsafe.Pointer(&A[0])),
        (*C.int)(unsafe.Pointer(&lda)),
        (unsafe.Pointer(&B[0])),
        (*C.int)(unsafe.Pointer(&ldb)),
        (unsafe.Pointer(&beta)),
        (unsafe.Pointer(&C[0])),
        (*C.int)(unsafe.Pointer(&ldc)))
}

func zhemm(side, uplo string, M int, N int,
    alpha complex128, A []complex128, lda int, B []complex128, ldb int, beta complex128,
    C []complex128, ldc int) {

    cuplo := C.CString(uplo)
    defer C.free(unsafe.Pointer(cuplo))
    cside := C.CString(side)
    defer C.free(unsafe.Pointer(cside))

    C.zhemm_(cside, cuplo,
        (*C.int)(unsafe.Pointer(&M)),
        (*C.int)(unsafe.Pointer(&N)),
        (unsafe.Pointer(&alpha)),
        (unsafe.Pointer(&A[0])),
        (*C.int)(unsafe.Pointer(&lda)),
        (unsafe.Pointer(&B[0])),
        (*C.int)(unsafe.Pointer(&ldb)),
        (unsafe.Pointer(&beta)),
        (unsafe.Pointer(&C[0])),
        (*C.int)(unsafe.Pointer(&ldc)))
}

func zsyrk(uplo, trans string, N int, K int,
    alpha complex128, A []complex128, lda int, beta complex128,
    C []complex128, ldc int) {

    cuplo := C.CString(uplo)
    defer C.free(unsafe.Pointer(cuplo))
    ctrans := C.CString(trans)
    defer C.free(unsafe.Pointer(ctrans))
    //cside := C.CString(side)
    //defer C.free(unsafe.Pointer(cside))

    C.zsyrk_(cuplo, ctrans,
        (*C.int)(unsafe.Pointer(&N)),
        (*C.int)(unsafe.Pointer(&K)),
        (unsafe.Pointer(&alpha)),
        (unsafe.Pointer(&A[0])),
        (*C.int)(unsafe.Pointer(&lda)),
        (unsafe.Pointer(&beta)),
        (unsafe.Pointer(&C[0])),
        (*C.int)(unsafe.Pointer(&ldc)))
}

func zherk(uplo, trans string, N int, K int,
    alpha complex128, A []complex128, lda int, beta float64,
    B []complex128, ldb int) {

    cuplo := C.CString(uplo)
    defer C.free(unsafe.Pointer(cuplo))
    ctrans := C.CString(trans)
    defer C.free(unsafe.Pointer(ctrans))
    //cside := C.CString(side)
    //defer C.free(unsafe.Pointer(cside))

    C.zherk_(cuplo, ctrans,
        (*C.int)(unsafe.Pointer(&N)),
        (*C.int)(unsafe.Pointer(&K)),
        (*C.double)(unsafe.Pointer(&alpha)),
        (unsafe.Pointer(&A[0])),
        (*C.int)(unsafe.Pointer(&lda)),
        (*C.double)(unsafe.Pointer(&beta)),
        (unsafe.Pointer(&B[0])),
        (*C.int)(unsafe.Pointer(&ldb)))
}

func zsyr2k(uplo, trans string, N int, K int,
    alpha complex128, A []complex128, lda int, B []complex128, ldb int, beta complex128,
    C []complex128, ldc int) {

    cuplo := C.CString(uplo)
    defer C.free(unsafe.Pointer(cuplo))
    ctrans := C.CString(trans)
    defer C.free(unsafe.Pointer(ctrans))
    //cside := C.CString(side)
    //defer C.free(unsafe.Pointer(cside))

    C.zsyr2k_(cuplo, ctrans,
        (*C.int)(unsafe.Pointer(&N)),
        (*C.int)(unsafe.Pointer(&K)),
        (unsafe.Pointer(&alpha)),
        (unsafe.Pointer(&A[0])),
        (*C.int)(unsafe.Pointer(&lda)),
        (unsafe.Pointer(&B[0])),
        (*C.int)(unsafe.Pointer(&ldb)),
        (unsafe.Pointer(&beta)),
        (unsafe.Pointer(&C[0])),
        (*C.int)(unsafe.Pointer(&ldc)))
}

func zher2k(uplo, trans string, N int, K int,
    alpha complex128, A []complex128, lda int, B []complex128, ldb int, beta float64,
    C []complex128, ldc int) {

    cuplo := C.CString(uplo)
    defer C.free(unsafe.Pointer(cuplo))
    ctrans := C.CString(trans)
    defer C.free(unsafe.Pointer(ctrans))
    //cside := C.CString(side)
    //defer C.free(unsafe.Pointer(cside))

    C.zsyr2k_(cuplo, ctrans,
        (*C.int)(unsafe.Pointer(&N)),
        (*C.int)(unsafe.Pointer(&K)),
        (unsafe.Pointer(&alpha)),
        (unsafe.Pointer(&A[0])),
        (*C.int)(unsafe.Pointer(&lda)),
        (unsafe.Pointer(&B[0])),
        (*C.int)(unsafe.Pointer(&ldb)),
        (unsafe.Pointer(&beta)),
        (unsafe.Pointer(&C[0])),
        (*C.int)(unsafe.Pointer(&ldc)))
}

func ztrmm(side, uplo, transA, diag string,
    M int, N int, alpha complex128, A []complex128, lda int, B []complex128, ldb int) {

    cuplo := C.CString(uplo)
    defer C.free(unsafe.Pointer(cuplo))
    ctrans := C.CString(transA)
    defer C.free(unsafe.Pointer(ctrans))
    cside := C.CString(side)
    defer C.free(unsafe.Pointer(cside))
    cdiag := C.CString(diag)
    defer C.free(unsafe.Pointer(cdiag))

    C.ztrmm_(cside, cuplo, ctrans, cdiag,
        (*C.int)(unsafe.Pointer(&M)),
        (*C.int)(unsafe.Pointer(&N)),
        (unsafe.Pointer(&alpha)),
        (unsafe.Pointer(&A[0])),
        (*C.int)(unsafe.Pointer(&lda)),
        (unsafe.Pointer(&B[0])),
        (*C.int)(unsafe.Pointer(&ldb)))
}

func ztrsm(side, uplo, transA, diag string,
    M int, N int, alpha complex128, A []complex128, lda int, B []complex128, ldb int) {

    cuplo := C.CString(uplo)
    defer C.free(unsafe.Pointer(cuplo))
    ctrans := C.CString(transA)
    defer C.free(unsafe.Pointer(ctrans))
    cside := C.CString(side)
    defer C.free(unsafe.Pointer(cside))
    cdiag := C.CString(diag)
    defer C.free(unsafe.Pointer(cdiag))

    C.ztrsm_(cside, cuplo, ctrans, cdiag,
        (*C.int)(unsafe.Pointer(&M)),
        (*C.int)(unsafe.Pointer(&N)),
        (unsafe.Pointer(&alpha)),
        (unsafe.Pointer(&A[0])),
        (*C.int)(unsafe.Pointer(&lda)),
        (unsafe.Pointer(&B[0])),
        (*C.int)(unsafe.Pointer(&ldb)))
}

// Local Variables:
// tab-width: 4
// End:
