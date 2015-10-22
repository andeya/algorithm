
// Copyright (c) Harri Rautila, 2012

// This file is part of github.com/hrautila/matrix package. It is free software,
// distributed under the terms of GNU Lesser General Public License Version 3, or
// any later version. See the COPYING tile included in this archive.

package calgo

/*

extern void matmult_vp_notrans(double *C, const double *A, const double *B,
			    const double alpha,
			    int M, int N, int P, int S, int L, int R, int E,
			    int vlen);

extern void matmult_vpur_notrans(double *C, const double *A, const double *B,
			    const double alpha, int ldC, int ldA, int ldB,
			    int M, int N, int P, int S, int L, int R, int E,
			    int vlen);

*/
// #cgo CFLAGS: -O3 -funroll-loops
import "C"
import "unsafe"

// Calculate C +=alpha*A*B where C is M*N, A is M*P, B is P*N and alpha is scalar.
// Arrays C, A, B are column major order matrix data arrays.
func Mult(C, A, B []float64, alpha float64, ldC, ldA, ldB, M, N, P int) {
    C.matmult_vpur_notrans(
        (*C.double)(unsafe.Pointer(&C[0])),
        (*C.double)(unsafe.Pointer(&A[0])),
        (*C.double)(unsafe.Pointer(&B[0])), C.double(alpha),
        C.int(ldC), C.int(ldA), C.int(ldB),
        C.int(M), C.int(N), C.int(P), C.int(0), C.int(N), C.int(0), C.int(M),
        C.int(0))
}


// C += alpha*A*B for block defined by rows [R:E] and columns [S:L].
// C is M*N, A is M*P, B is P*N and 0 <= L < S < N and 0 <= R < E < M
func BlkMult(C, A, B []float64, alpha float64, ldC, ldA, ldB, M, N, P, S, L, R, E int) {
    C.matmult_vpur_notrans(
        (*C.double)(unsafe.Pointer(&C[0])),
        (*C.double)(unsafe.Pointer(&A[0])),
        (*C.double)(unsafe.Pointer(&B[0])), C.double(alpha),
        C.int(ldC), C.int(ldA), C.int(ldB),
        C.int(M), C.int(N), C.int(P), C.int(S), C.int(L), C.int(R), C.int(E),
        C.int(0))
}

// C += alpha*A*B for block defined by rows [R:E] and columns [S:L].
// C is M*N, A is M*P, B is P*N and 0 <= L < S < N and 0 <= R < E < M
// Panels in A and B are accessed in vlen blocks and accumulated to C.
func BlkMultVp(C, A, B []float64, alpha float64, ldC, ldA, ldB, M, N, P, S, L, R, E, vlen int) {
    C.matmult_vpur_notrans(
        (*C.double)(unsafe.Pointer(&C[0])),
        (*C.double)(unsafe.Pointer(&A[0])),
        (*C.double)(unsafe.Pointer(&B[0])),
        C.double(alpha), C.int(ldC), C.int(ldA), C.int(ldB),
        C.int(M), C.int(N), C.int(P), C.int(S), C.int(L), C.int(R), C.int(E),
        C.int(vlen))
}


// Local Variables:
// tab-width: 4
// indent-tabs-mode: nil
// End:

