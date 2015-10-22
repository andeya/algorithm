// Copyright (c) Harri Rautila, 2012,2013

// This file is part of github.com/hrautila/matops package. It is free software,
// distributed under the terms of GNU Lesser General Public License Version 3, or
// any later version. See the COPYING tile included in this archive.

package calgo

// -O3 -msse4.1 -funroll-loops -fomit-frame-pointer -ffast-math

// #cgo CFLAGS: -O3 -msse4.1 -fomit-frame-pointer -ffast-math
// #cgo LDFLAGS: -lm
// #include "cmops.h"
import "C"
import "unsafe"

type Flags int

const (
	TRANSA  = (1 << iota) // 0x1  ; transpose A
	TRANSB                // 0x2  ; transpose B
	LOWER                 // 0x4  ; lower tridiagonal
	UPPER                 // 0x8  ; upper tridiagonal
	LEFT                  // 0x10 ; A on left side
	RIGHT                 // 0x20 ; A on right side
	UNIT                  // 0x40 ; unit diagonal
	TRANS                 // 0x80 ; generic transpose
	NOTRANS = 0
	NULL    = 0
)

// Matrix-Matrix operators

// matrix-matrix: A = alpha*A + beta*B
func DScalePlus(A, B []float64, alpha, beta float64, flags Flags, ldA, ldB, S, L, R, E int) {
	var Am C.mdata_t
	var Bm C.mdata_t

	if B == nil || A == nil {
		return
	}
	Am.md = (*C.double)(unsafe.Pointer(&A[0]))
	Am.step = C.int(ldA)
	Bm.md = (*C.double)(unsafe.Pointer(&B[0]))
	Bm.step = C.int(ldB)

	C.dmmat_scale_plus(
		(*C.mdata_t)(unsafe.Pointer(&Am)),
		(*C.mdata_t)(unsafe.Pointer(&Bm)),
		C.double(alpha), C.double(beta), C.int(flags),
		C.int(S), C.int(L),
		C.int(R), C.int(E))
}

// Generic matrix-matrix multiplication for block [R:E, S:L] with panel length P.
//
// if trans is NOTRANS then calculates
//   C = alpha*A*B + beta*C; C is M*N, A is M*P and B is P*N;
// if trans is TRANSA then calculates
//   C = alpha*A.T*B + beta*C; C is M*N, A is P*M and B is P*N;
// if trans is TRANSB then calculates
//   C = alpha*A*B.T + beta*C; C is M*N, A is M*P and B is N*P;
// if trans is TRANSA|TRANSB then calculates
//   C = alpha*A.T*B.T + beta*C; C is M*N, A is P*M and B is N*P;
//
func DMult(C, A, B []float64, alpha, beta float64, trans Flags, ldC, ldA, ldB, P, S, L, R, E, H, NB, MB int) {

	var Cm C.mdata_t
	var Am C.mdata_t
	var Bm C.mdata_t

	if C == nil || B == nil || A == nil {
		return
	}
	Cm.md = (*C.double)(unsafe.Pointer(&C[0]))
	Cm.step = C.int(ldC)
	Am.md = (*C.double)(unsafe.Pointer(&A[0]))
	Am.step = C.int(ldA)
	Bm.md = (*C.double)(unsafe.Pointer(&B[0]))
	Bm.step = C.int(ldB)

	C.dmult_mm_blocked3(
		(*C.mdata_t)(unsafe.Pointer(&Cm)),
		(*C.mdata_t)(unsafe.Pointer(&Am)),
		(*C.mdata_t)(unsafe.Pointer(&Bm)),
		C.double(alpha), C.double(beta), C.int(trans),
		C.int(P), C.int(S), C.int(L), C.int(R), C.int(E),
		C.int(H), C.int(NB), C.int(MB))
}

//
func DMultSymm(C, A, B []float64, alpha, beta float64, flags Flags, ldC, ldA, ldB, P, S, L, R, E, H, NB, MB int) {

	var Cm C.mdata_t
	var Am C.mdata_t
	var Bm C.mdata_t

	if C == nil || B == nil || A == nil {
		return
	}
	Cm.md = (*C.double)(unsafe.Pointer(&C[0]))
	Cm.step = C.int(ldC)
	Am.md = (*C.double)(unsafe.Pointer(&A[0]))
	Am.step = C.int(ldA)
	Bm.md = (*C.double)(unsafe.Pointer(&B[0]))
	Bm.step = C.int(ldB)

	C.dmult_symm_blocked2(
		(*C.mdata_t)(unsafe.Pointer(&Cm)),
		(*C.mdata_t)(unsafe.Pointer(&Am)),
		(*C.mdata_t)(unsafe.Pointer(&Bm)),
		C.double(alpha), C.double(beta), C.int(flags),
		C.int(P), C.int(S), C.int(L), C.int(R), C.int(E),
		C.int(H), C.int(NB), C.int(MB))
}

// blas TRMM; unblocked
// S is the start column (LEFT), row (RIGHT); E is the end column (LEFT), row (RIGHT)
func DTrmmUnblk(B, A []float64, alpha float64, flags Flags, ldB, ldA, N, S, E, NB int) {
	var Bm C.mdata_t
	var Am C.mdata_t

	if B == nil || A == nil {
		return
	}
	if N == 0 || E-S <= 0 {
		return
	}
	Bm.md = (*C.double)(unsafe.Pointer(&B[0]))
	Bm.step = C.int(ldB)
	Am.md = (*C.double)(unsafe.Pointer(&A[0]))
	Am.step = C.int(ldA)

	C.dmmat_trid_unb(
		(*C.mdata_t)(unsafe.Pointer(&Bm)),
		(*C.mdata_t)(unsafe.Pointer(&Am)),
		C.double(alpha), C.int(flags), C.int(N), C.int(S), C.int(E),
		/*C.int(R), C.int(E),*/)

}

// blas TRMM; blocked
func DTrmmBlk(B, A []float64, alpha float64, flags Flags, ldB, ldA, N, S, E, NB int) {
	var Bm C.mdata_t
	var Am C.mdata_t

	if B == nil || A == nil {
		return
	}
	if N == 0 || E-S <= 0 {
		return
	}
	Bm.md = (*C.double)(unsafe.Pointer(&B[0]))
	Bm.step = C.int(ldB)
	Am.md = (*C.double)(unsafe.Pointer(&A[0]))
	Am.step = C.int(ldA)

	C.dmmat_trmm_blk(
		(*C.mdata_t)(unsafe.Pointer(&Bm)),
		(*C.mdata_t)(unsafe.Pointer(&Am)),
		C.double(alpha), C.int(flags), C.int(N), C.int(S), C.int(E),
		/*C.int(R), C.int(E),*/ C.int(NB))

}

// blas TRSM; unblocked
// S is the start column (LEFT), row (RIGHT); E is the end column (LEFT), row (RIGHT)
func DSolveUnblk(B, A []float64, alpha float64, flags Flags, ldB, ldA, N, S, E int) {
	var Bm C.mdata_t
	var Am C.mdata_t

	if B == nil || A == nil {
		return
	}
	if N == 0 || E-S <= 0 {
		return
	}
	Bm.md = (*C.double)(unsafe.Pointer(&B[0]))
	Bm.step = C.int(ldB)
	Am.md = (*C.double)(unsafe.Pointer(&A[0]))
	Am.step = C.int(ldA)

	C.dmmat_solve_unb(
		(*C.mdata_t)(unsafe.Pointer(&Bm)),
		(*C.mdata_t)(unsafe.Pointer(&Am)),
		C.double(alpha), C.int(flags), C.int(N), C.int(S), C.int(E))

}

// blas TRSM; blocked
// S is the start column (LEFT), row (RIGHT); E is the end column (LEFT), row (RIGHT)
func DSolveBlk(B, A []float64, alpha float64, flags Flags, ldB, ldA, N, S, E, NB int) {
	var Bm C.mdata_t
	var Am C.mdata_t

	if B == nil || A == nil {
		return
	}
	if N == 0 || E-S <= 0 {
		return
	}
	Bm.md = (*C.double)(unsafe.Pointer(&B[0]))
	Bm.step = C.int(ldB)
	Am.md = (*C.double)(unsafe.Pointer(&A[0]))
	Am.step = C.int(ldA)

	C.dmmat_solve_blk(
		(*C.mdata_t)(unsafe.Pointer(&Bm)),
		(*C.mdata_t)(unsafe.Pointer(&Am)),
		C.double(alpha), C.int(flags), C.int(N), C.int(S), C.int(E), C.int(NB))

}

// blas SYRK; blocked
// S is the start column and row in C; E is the end column and row in C
func DSymmRankBlk(C, A []float64, alpha, beta float64, flags Flags, ldC, ldA, N, S, E, H, NB int) {
	var Cm C.mdata_t
	var Am C.mdata_t

	if C == nil || A == nil {
		return
	}
	if N == 0 || E-S <= 0 {
		return
	}
	Cm.md = (*C.double)(unsafe.Pointer(&C[0]))
	Cm.step = C.int(ldC)
	Am.md = (*C.double)(unsafe.Pointer(&A[0]))
	Am.step = C.int(ldA)

	C.dmmat_rank_blk(
		(*C.mdata_t)(unsafe.Pointer(&Cm)),
		(*C.mdata_t)(unsafe.Pointer(&Am)),
		C.double(alpha), C.double(beta),
		C.int(flags), C.int(N), C.int(S), C.int(E), C.int(H), C.int(NB))

}

// blas SYR2K; blocked
// S is the start column and row in C; E is the end column and row in C
func DSymmRank2Blk(C, A, B []float64, alpha, beta float64, flags Flags, ldC, ldA, ldB, N, S, E, H, NB int) {
	var Cm C.mdata_t
	var Am C.mdata_t
	var Bm C.mdata_t

	if C == nil || A == nil || B == nil {
		return
	}
	if N == 0 || E-S <= 0 {
		return
	}
	Cm.md = (*C.double)(unsafe.Pointer(&C[0]))
	Cm.step = C.int(ldC)
	Am.md = (*C.double)(unsafe.Pointer(&A[0]))
	Am.step = C.int(ldA)
	Bm.md = (*C.double)(unsafe.Pointer(&B[0]))
	Bm.step = C.int(ldB)

	C.dmmat_rank2_blk(
		(*C.mdata_t)(unsafe.Pointer(&Cm)),
		(*C.mdata_t)(unsafe.Pointer(&Am)),
		(*C.mdata_t)(unsafe.Pointer(&Bm)),
		C.double(alpha), C.double(beta),
		C.int(flags), C.int(N), C.int(S), C.int(E), C.int(H), C.int(NB))

}

// Generic triangular matrix update; blocked
// S is the start column and row in C; E is the end column and row in C
func DTrmUpdBlk(C, A, B []float64, alpha, beta float64, flags Flags, ldC, ldA, ldB, N, S, E, H, NB int) {
	var Cm C.mdata_t
	var Am C.mdata_t
	var Bm C.mdata_t

	if C == nil || A == nil || B == nil {
		return
	}
	if N == 0 || E-S <= 0 {
		return
	}
	Cm.md = (*C.double)(unsafe.Pointer(&C[0]))
	Cm.step = C.int(ldC)
	Am.md = (*C.double)(unsafe.Pointer(&A[0]))
	Am.step = C.int(ldA)
	Bm.md = (*C.double)(unsafe.Pointer(&B[0]))
	Bm.step = C.int(ldB)

	C.dmmat_trmupd_blk(
		(*C.mdata_t)(unsafe.Pointer(&Cm)),
		(*C.mdata_t)(unsafe.Pointer(&Am)),
		(*C.mdata_t)(unsafe.Pointer(&Bm)),
		C.double(alpha), C.double(beta),
		C.int(flags), C.int(N), C.int(S), C.int(E), C.int(H), C.int(NB))

}

// Matrix-Vector operators

// blas GEMV; blocked version
// Y = alpha*A*X + beta*Y; Y is M*1, X is N*1 and A is M*N
func DMultMV(Y, A, X []float64, alpha, beta float64, flags Flags, incY, ldA, incX, S, L, R, E, H, MB int) {
	var Yv C.mvec_t
	var Xv C.mvec_t
	var Am C.mdata_t

	if Y == nil || A == nil || X == nil {
		return
	}
	Yv.md = (*C.double)(unsafe.Pointer(&Y[0]))
	Yv.inc = C.int(incY)
	Xv.md = (*C.double)(unsafe.Pointer(&X[0]))
	Xv.inc = C.int(incX)
	Am.md = (*C.double)(unsafe.Pointer(&A[0]))
	Am.step = C.int(ldA)

	C.dmult_gemv_blocked(
		(*C.mvec_t)(unsafe.Pointer(&Yv)),
		(*C.mdata_t)(unsafe.Pointer(&Am)),
		(*C.mvec_t)(unsafe.Pointer(&Xv)),
		C.double(alpha), C.double(beta), C.int(flags),
		C.int(S), C.int(L), C.int(R), C.int(E),
		C.int(H), C.int(MB))
}

// blas GER; blocked version
// A = A + alpha * x * y.T; A is M*N, x is M*1, Y is N*1, 0 < R < E <= M, 0 < S < L <= N
func DRankMV(A, X, Y []float64, alpha float64, ldA, incX, incY, S, L, R, E, NB, MB int) {
	var Yv C.mvec_t
	var Xv C.mvec_t
	var Am C.mdata_t

	if A == nil || X == nil || Y == nil {
		return
	}

	Yv.md = (*C.double)(unsafe.Pointer(&Y[0]))
	Yv.inc = C.int(incY)
	Xv.md = (*C.double)(unsafe.Pointer(&X[0]))
	Xv.inc = C.int(incX)
	Am.md = (*C.double)(unsafe.Pointer(&A[0]))
	Am.step = C.int(ldA)

	C.dmvec_rank(
		(*C.mdata_t)(unsafe.Pointer(&Am)),
		(*C.mvec_t)(unsafe.Pointer(&Xv)),
		(*C.mvec_t)(unsafe.Pointer(&Yv)),
		C.double(alpha),
		C.int(S), C.int(L), C.int(R), C.int(E),
		C.int(NB), C.int(MB))
}

// A = A + alpha * x * y.T; A is M*N, x is M*1, Y is N*1, 0 < R < E <= M, 0 < S < L <= N

// blas SYR; blocked version
func DSymmRankMV(A, X []float64, alpha float64, flags Flags, ldA, incX, S, L, NB int) {
	var Xv C.mvec_t
	var Am C.mdata_t

	if A == nil || X == nil {
		return
	}
	Xv.md = (*C.double)(unsafe.Pointer(&X[0]))
	Xv.inc = C.int(incX)
	Am.md = (*C.double)(unsafe.Pointer(&A[0]))
	Am.step = C.int(ldA)

	C.dmvec_symv_rank(
		(*C.mdata_t)(unsafe.Pointer(&Am)),
		(*C.mvec_t)(unsafe.Pointer(&Xv)),
		C.double(alpha), C.int(flags),
		C.int(S), C.int(L), C.int(NB))
}

// A = A + alpha * x * y.T; A is M*N, x is M*1, Y is N*1, 0 < R < E <= M, 0 < S < L <= N

// blas SYR2; blocked version
func DSymmRank2MV(A, X, Y []float64, alpha float64, flags Flags, ldA, incX, incY, S, L, NB int) {
	var Xv C.mvec_t
	var Yv C.mvec_t
	var Am C.mdata_t

	if A == nil || X == nil || Y == nil {
		return
	}
	Xv.md = (*C.double)(unsafe.Pointer(&X[0]))
	Xv.inc = C.int(incX)
	Yv.md = (*C.double)(unsafe.Pointer(&Y[0]))
	Yv.inc = C.int(incY)
	Am.md = (*C.double)(unsafe.Pointer(&A[0]))
	Am.step = C.int(ldA)

	C.dmvec_symv_rank2(
		(*C.mdata_t)(unsafe.Pointer(&Am)),
		(*C.mvec_t)(unsafe.Pointer(&Xv)),
		(*C.mvec_t)(unsafe.Pointer(&Yv)),
		C.double(alpha), C.int(flags),
		C.int(S), C.int(L), C.int(NB))
}

// generic triangular matrix rank update; A = A + alpha*X*Y.T
func DTrmUpdMV(A, X, Y []float64, alpha float64, flags Flags, ldA, incX, incY, S, L, NB int) {
	var Xv C.mvec_t
	var Yv C.mvec_t
	var Am C.mdata_t

	if A == nil || X == nil || Y == nil {
		return
	}
	Xv.md = (*C.double)(unsafe.Pointer(&X[0]))
	Xv.inc = C.int(incX)
	Yv.md = (*C.double)(unsafe.Pointer(&Y[0]))
	Yv.inc = C.int(incY)
	Am.md = (*C.double)(unsafe.Pointer(&A[0]))
	Am.step = C.int(ldA)

	C.dmvec_trmv_upd(
		(*C.mdata_t)(unsafe.Pointer(&Am)),
		(*C.mvec_t)(unsafe.Pointer(&Xv)),
		(*C.mvec_t)(unsafe.Pointer(&Yv)),
		C.double(alpha), C.int(flags),
		C.int(S), C.int(L), C.int(NB))
}

// blas TSMV; unblocked version
func DSolveUnblkMV(X, A []float64, flags Flags, incX, ldA, N int) {
	var Xv C.mvec_t
	var Am C.mdata_t

	if A == nil || X == nil {
		return
	}
	Xv.md = (*C.double)(unsafe.Pointer(&X[0]))
	Xv.inc = C.int(incX)
	Am.md = (*C.double)(unsafe.Pointer(&A[0]))
	Am.step = C.int(ldA)

	C.dmvec_solve_unb(
		(*C.mvec_t)(unsafe.Pointer(&Xv)),
		(*C.mdata_t)(unsafe.Pointer(&Am)),
		C.int(flags), C.int(N))

}

// blas TSMV; blocked version
func DSolveBlkMV(X, A []float64, flags Flags, incX, ldA, N, NB int) {
	var Xv C.mvec_t
	var Am C.mdata_t

	if A == nil || X == nil {
		return
	}
	Xv.md = (*C.double)(unsafe.Pointer(&X[0]))
	Xv.inc = C.int(incX)
	Am.md = (*C.double)(unsafe.Pointer(&A[0]))
	Am.step = C.int(ldA)

	C.dmvec_solve_blocked(
		(*C.mvec_t)(unsafe.Pointer(&Xv)),
		(*C.mdata_t)(unsafe.Pointer(&Am)),
		C.int(flags), C.int(N), C.int(NB))

}

// blas TRMV; unblocked
func DTrimvUnblkMV(X, A []float64, flags Flags, incX, ldA, N int) {
	var Xv C.mvec_t
	var Am C.mdata_t

	if A == nil || X == nil {
		return
	}
	Xv.md = (*C.double)(unsafe.Pointer(&X[0]))
	Xv.inc = C.int(incX)
	Am.md = (*C.double)(unsafe.Pointer(&A[0]))
	Am.step = C.int(ldA)

	C.dmvec_trid_unb(
		(*C.mvec_t)(unsafe.Pointer(&Xv)),
		(*C.mdata_t)(unsafe.Pointer(&Am)),
		C.int(flags), C.int(N))

}

// Z[0] = beta*Z[0] + alpha * X * Y
func DDotSum(Z, X, Y []float64, alpha, beta float64, incZ, incX, incY, N int) {

	var Zv C.mvec_t
	var Xv C.mvec_t
	var Yv C.mvec_t

	if Z == nil || X == nil || Y == nil || N <= 0 {
		return
	}
	Xv.md = (*C.double)(unsafe.Pointer(&X[0]))
	Xv.inc = C.int(incX)
	Zv.md = (*C.double)(unsafe.Pointer(&Z[0]))
	Zv.inc = C.int(incZ)
	Yv.md = (*C.double)(unsafe.Pointer(&Y[0]))
	Yv.inc = C.int(incY)

	C.dvec_dots(
		(*C.mvec_t)(unsafe.Pointer(&Zv)),
		(*C.mvec_t)(unsafe.Pointer(&Xv)),
		(*C.mvec_t)(unsafe.Pointer(&Yv)),
		C.double(alpha), C.double(beta),
		C.int(N))
}

// return: alpha * X * Y
func DDot(X, Y []float64, alpha float64, incX, incY, N int) float64 {

	var dot C.double
	var Xv C.mvec_t
	var Yv C.mvec_t

	if X == nil || Y == nil || N <= 0 {
		return 0.0
	}
	Xv.md = (*C.double)(unsafe.Pointer(&X[0]))
	Xv.inc = C.int(incX)
	Yv.md = (*C.double)(unsafe.Pointer(&Y[0]))
	Yv.inc = C.int(incY)

	dot = C.dvec_dot(
		(*C.mvec_t)(unsafe.Pointer(&Xv)),
		(*C.mvec_t)(unsafe.Pointer(&Yv)),
		C.double(alpha),
		C.int(N))
	return float64(dot)
}

// Y := alpha*X + Y
func DAxpy(Y, X []float64, alpha float64, incX, incY, N int) {

	var Xv C.mvec_t
	var Yv C.mvec_t

	if X == nil || Y == nil || N <= 0 {
		return
	}
	Xv.md = (*C.double)(unsafe.Pointer(&X[0]))
	Xv.inc = C.int(incX)
	Yv.md = (*C.double)(unsafe.Pointer(&Y[0]))
	Yv.inc = C.int(incY)

	C.dvec_axpy(
		(*C.mvec_t)(unsafe.Pointer(&Yv)),
		(*C.mvec_t)(unsafe.Pointer(&Xv)),
		C.double(alpha),
		C.int(N))
}

// return: sum (abs(X[i]-Y[i]))^2
func DiffNorm2(X, Y []float64, incX, incY, N int) float64 {

	var nrm C.double
	var Xv C.mvec_t
	var Yv C.mvec_t

	if X == nil || Y == nil || N <= 0 {
		return 0.0
	}
	Xv.md = (*C.double)(unsafe.Pointer(&X[0]))
	Xv.inc = C.int(incX)
	Yv.md = (*C.double)(unsafe.Pointer(&Y[0]))
	Yv.inc = C.int(incY)

	nrm = C.dvec_diff_nrm2(
		(*C.mvec_t)(unsafe.Pointer(&Xv)),
		(*C.mvec_t)(unsafe.Pointer(&Yv)),
		C.int(N))
	return float64(nrm)
}

// return: sum (abs(X[i]))^2; Euclidaen norm
func DNorm2(X []float64, incX, N int) float64 {

	var nrm C.double
	var Xv C.mvec_t

	if X == nil || N <= 0 {
		return 0.0
	}
	Xv.md = (*C.double)(unsafe.Pointer(&X[0]))
	Xv.inc = C.int(incX)

	nrm = C.dvec_nrm2(
		(*C.mvec_t)(unsafe.Pointer(&Xv)),
		C.int(N))
	return float64(nrm)
}

// return: sum (abs(X[i]))
func DAsum(X []float64, incX, N int) float64 {

	var asum C.double
	var Xv C.mvec_t

	if X == nil || N <= 0 {
		return 0.0
	}
	Xv.md = (*C.double)(unsafe.Pointer(&X[0]))
	Xv.inc = C.int(incX)

	asum = C.dvec_asum(
		(*C.mvec_t)(unsafe.Pointer(&Xv)),
		C.int(N))
	return float64(asum)
}

// return: index of max absolute value
func DIAMax(X []float64, incX, N int) int {

	var ix C.int
	var Xv C.mvec_t

	if X == nil || N <= 0 {
		return -1
	}
	Xv.md = (*C.double)(unsafe.Pointer(&X[0]))
	Xv.inc = C.int(incX)

	ix = C.dvec_iamax(
		(*C.mvec_t)(unsafe.Pointer(&Xv)),
		C.int(N))
	return int(ix)
}

func DSwap(X, Y []float64, incX, incY, N int) {

	var Xv C.mvec_t
	var Yv C.mvec_t

	if X == nil || Y == nil || N <= 0 {
		return
	}
	Xv.md = (*C.double)(unsafe.Pointer(&X[0]))
	Xv.inc = C.int(incX)
	Yv.md = (*C.double)(unsafe.Pointer(&Y[0]))
	Yv.inc = C.int(incY)

	C.dvec_swap(
		(*C.mvec_t)(unsafe.Pointer(&Xv)),
		(*C.mvec_t)(unsafe.Pointer(&Yv)),
		C.int(N))
}

// copying: X := Y
func DCopy(X, Y []float64, incX, incY, N int) {

	var Xv C.mvec_t
	var Yv C.mvec_t

	if X == nil || Y == nil || N <= 0 {
		return
	}
	Xv.md = (*C.double)(unsafe.Pointer(&X[0]))
	Xv.inc = C.int(incX)
	Yv.md = (*C.double)(unsafe.Pointer(&Y[0]))
	Yv.inc = C.int(incY)

	C.dvec_copy(
		(*C.mvec_t)(unsafe.Pointer(&Xv)),
		(*C.mvec_t)(unsafe.Pointer(&Yv)),
		C.int(N))
}

// inverse scaling: X = X/alpha
func DInvScal(X []float64, alpha float64, incX, N int) {

	var Xv C.mvec_t

	if X == nil || N <= 0 {
		return
	}
	Xv.md = (*C.double)(unsafe.Pointer(&X[0]))
	Xv.inc = C.int(incX)

	C.dvec_invscal(
		(*C.mvec_t)(unsafe.Pointer(&Xv)),
		C.double(alpha), C.int(N))
}

// scaling: X = alpha*X
func DScal(X []float64, alpha float64, incX, N int) {

	var Xv C.mvec_t

	if X == nil || N <= 0 {
		return
	}
	Xv.md = (*C.double)(unsafe.Pointer(&X[0]))
	Xv.inc = C.int(incX)

	C.dvec_scal(
		(*C.mvec_t)(unsafe.Pointer(&Xv)),
		C.double(alpha), C.int(N))
}

// Local Variables:
// tab-width: 4
// indent-tabs-mode: nil
// End:
