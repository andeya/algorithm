

package calgo

// #cgo CFLAGS: -O3 -msse4.1 -funroll-loops -fomit-frame-pointer -ffast-math 
// #include "cmops.h"
import "C"
import "unsafe"

func DMult0(C, A, B []float64, alpha, beta float64, trans Flags, ldC, ldA, ldB, P, S, L, R, E, H, NB, MB int) {

    var Cm C.mdata_t
    var Am C.mdata_t
    var Bm C.mdata_t

    if C == nil || B == nil || A == nil {
        return
    }
    Cm.md =  (*C.double)(unsafe.Pointer(&C[0]))
    Cm.step = C.int(ldC)
    Am.md =  (*C.double)(unsafe.Pointer(&A[0]))
    Am.step = C.int(ldA)
    Bm.md =  (*C.double)(unsafe.Pointer(&B[0]))
    Bm.step = C.int(ldB)

    C.dmult_mm_blocked2(
        (*C.mdata_t)(unsafe.Pointer(&Cm)),
        (*C.mdata_t)(unsafe.Pointer(&Am)),
        (*C.mdata_t)(unsafe.Pointer(&Bm)),
        C.double(alpha), C.double(beta), C.int(trans),
        C.int(P), C.int(S), C.int(L), C.int(R), C.int(E),
        C.int(H), C.int(NB), C.int(MB))
}

func DMultOld(C, A, B []float64, alpha, beta float64, trans Flags, ldC, ldA, ldB, P, S, L, R, E, H, NB, MB int) {

    var Cm C.mdata_t
    var Am C.mdata_t
    var Bm C.mdata_t

    Cm.md =  (*C.double)(unsafe.Pointer(&C[0]))
    Cm.step = C.int(ldC)
    Am.md =  (*C.double)(unsafe.Pointer(&A[0]))
    Am.step = C.int(ldA)
    Bm.md =  (*C.double)(unsafe.Pointer(&B[0]))
    Bm.step = C.int(ldB)

    C.dmult_mm_blocked(
        (*C.mdata_t)(unsafe.Pointer(&Cm)),
        (*C.mdata_t)(unsafe.Pointer(&Am)),
        (*C.mdata_t)(unsafe.Pointer(&Bm)),
        C.double(alpha), C.double(beta), C.int(trans),
        C.int(P), C.int(S), C.int(L), C.int(R), C.int(E),
        C.int(H), C.int(NB), C.int(MB))
}

func DMultSymmOld(C, A, B []float64, alpha, beta float64, flags Flags, ldC, ldA, ldB, P, S, L, R, E, H, NB, MB int) {

    var Cm C.mdata_t
    var Am C.mdata_t
    var Bm C.mdata_t

    Cm.md =  (*C.double)(unsafe.Pointer(&C[0]))
    Cm.step = C.int(ldC)
    Am.md =  (*C.double)(unsafe.Pointer(&A[0]))
    Am.step = C.int(ldA)
    Bm.md =  (*C.double)(unsafe.Pointer(&B[0]))
    Bm.step = C.int(ldB)

    C.dmult_symm_blocked(
        (*C.mdata_t)(unsafe.Pointer(&Cm)),
        (*C.mdata_t)(unsafe.Pointer(&Am)),
        (*C.mdata_t)(unsafe.Pointer(&Bm)),
        C.double(alpha), C.double(beta), C.int(flags),
        C.int(P), C.int(S), C.int(L), C.int(R), C.int(E),
        C.int(H), C.int(NB), C.int(MB))
}

func DSolveLower(X, A []float64, unit bool, incX, ldA, N, NB int) {
    var Xv C.mvec_t
    var Am C.mdata_t
    Xv.md =  (*C.double)(unsafe.Pointer(&X[0]))
    Xv.inc = C.int(incX)
    Am.md =  (*C.double)(unsafe.Pointer(&A[0]))
    Am.step = C.int(ldA)

    var flags Flags = LOWER
    if unit {
        flags |= UNIT
    }
    C.dmvec_solve_unb(
        (*C.mvec_t)(unsafe.Pointer(&Xv)),
        (*C.mdata_t)(unsafe.Pointer(&Am)),
        C.int(flags), C.int(N))

}

func DSolveLowerBlocked(X, A []float64, unit bool, incX, ldA, N, NB int) {
    var Xv C.mvec_t
    var Am C.mdata_t
    Xv.md =  (*C.double)(unsafe.Pointer(&X[0]))
    Xv.inc = C.int(incX)
    Am.md =  (*C.double)(unsafe.Pointer(&A[0]))
    Am.step = C.int(ldA)

    var flags Flags = LOWER
    if unit {
        flags |= UNIT
    }
    C.dmvec_solve_blocked(
        (*C.mvec_t)(unsafe.Pointer(&Xv)),
        (*C.mdata_t)(unsafe.Pointer(&Am)),
        C.int(flags), C.int(N), C.int(NB))

}

func DSolveUpper(X, A []float64, unit bool, incX, ldA, N, NB int) {
    var Xv C.mvec_t
    var Am C.mdata_t
    Xv.md =  (*C.double)(unsafe.Pointer(&X[0]))
    Xv.inc = C.int(incX)
    Am.md =  (*C.double)(unsafe.Pointer(&A[0]))
    Am.step = C.int(ldA)

    var flags Flags = UPPER
    if unit {
        flags |= UNIT
    }
    C.dmvec_solve_unb(
        (*C.mvec_t)(unsafe.Pointer(&Xv)),
        (*C.mdata_t)(unsafe.Pointer(&Am)),
        C.int(flags), C.int(N))

}

func DSolveUpperBlocked(X, A []float64, unit bool, incX, ldA, N, NB int) {
    var Xv C.mvec_t
    var Am C.mdata_t
    Xv.md =  (*C.double)(unsafe.Pointer(&X[0]))
    Xv.inc = C.int(incX)
    Am.md =  (*C.double)(unsafe.Pointer(&A[0]))
    Am.step = C.int(ldA)

    var flags Flags = UPPER
    if unit {
        flags |= UNIT
    }
    C.dmvec_solve_blocked(
        (*C.mvec_t)(unsafe.Pointer(&Xv)),
        (*C.mdata_t)(unsafe.Pointer(&Am)),
        C.int(flags), C.int(N), C.int(NB))

}

func DTrimvUpper(X, A []float64, unit bool, incX, ldA, N, NB int) {
    var Xv C.mvec_t
    var Am C.mdata_t
    Xv.md =  (*C.double)(unsafe.Pointer(&X[0]))
    Xv.inc = C.int(incX)
    Am.md =  (*C.double)(unsafe.Pointer(&A[0]))
    Am.step = C.int(ldA)

    var flags Flags = UPPER
    if unit {
        flags |= UNIT
    }
    C.dmvec_trid_unb(
        (*C.mvec_t)(unsafe.Pointer(&Xv)),
        (*C.mdata_t)(unsafe.Pointer(&Am)),
        C.int(flags), C.int(N))

}

func DTrimvUpperTransA(X, A []float64, unit bool, incX, ldA, N, NB int) {
    var Xv C.mvec_t
    var Am C.mdata_t
    Xv.md =  (*C.double)(unsafe.Pointer(&X[0]))
    Xv.inc = C.int(incX)
    Am.md =  (*C.double)(unsafe.Pointer(&A[0]))
    Am.step = C.int(ldA)

    var flags Flags = UPPER|TRANSA
    if unit {
        flags |= UNIT
    }
    C.dmvec_trid_unb(
        (*C.mvec_t)(unsafe.Pointer(&Xv)),
        (*C.mdata_t)(unsafe.Pointer(&Am)),
        C.int(flags), C.int(N))

}

func DTrimvLower(X, A []float64, unit bool, incX, ldA, N, NB int) {
    var Xv C.mvec_t
    var Am C.mdata_t
    Xv.md =  (*C.double)(unsafe.Pointer(&X[0]))
    Xv.inc = C.int(incX)
    Am.md =  (*C.double)(unsafe.Pointer(&A[0]))
    Am.step = C.int(ldA)

    var flags Flags = LOWER
    if unit {
        flags |= UNIT
    }
    C.dmvec_trid_unb(
        (*C.mvec_t)(unsafe.Pointer(&Xv)),
        (*C.mdata_t)(unsafe.Pointer(&Am)), 
        C.int(flags), C.int(N))

}

func DTrimvLowerTransA(X, A []float64, unit bool, incX, ldA, N, NB int) {
    var Xv C.mvec_t
    var Am C.mdata_t
    Xv.md =  (*C.double)(unsafe.Pointer(&X[0]))
    Xv.inc = C.int(incX)
    Am.md =  (*C.double)(unsafe.Pointer(&A[0]))
    Am.step = C.int(ldA)

    var flags Flags = LOWER|TRANSA
    if unit {
        flags |= UNIT
    }
    C.dmvec_trid_unb(
        (*C.mvec_t)(unsafe.Pointer(&Xv)),
        (*C.mdata_t)(unsafe.Pointer(&Am)), 
        C.int(flags), C.int(N))

}

func DTrmmUpper(B, A []float64, alpha float64, unit bool, ldB, ldA, N, S, L int) {
    var Bm C.mdata_t
    var Am C.mdata_t
    Bm.md =  (*C.double)(unsafe.Pointer(&B[0]))
    Bm.step = C.int(ldB)
    Am.md =  (*C.double)(unsafe.Pointer(&A[0]))
    Am.step = C.int(ldA)

    var flags Flags = UPPER
    if unit {
        flags |= UNIT
    }
    C.dmmat_trid_unb(
        (*C.mdata_t)(unsafe.Pointer(&Bm)),
        (*C.mdata_t)(unsafe.Pointer(&Am)),
        C.double(alpha), C.int(flags), C.int(N), C.int(S), C.int(L))

}

func DTrmmUpperTransA(B, A []float64, alpha float64, unit bool, ldB, ldA, N, S, L int) {
    var Bm C.mdata_t
    var Am C.mdata_t
    Bm.md =  (*C.double)(unsafe.Pointer(&B[0]))
    Bm.step = C.int(ldB)
    Am.md =  (*C.double)(unsafe.Pointer(&A[0]))
    Am.step = C.int(ldA)

    var flags Flags = UPPER|TRANSA
    if unit {
        flags |= UNIT
    }
    C.dmmat_trid_unb(
        (*C.mdata_t)(unsafe.Pointer(&Bm)),
        (*C.mdata_t)(unsafe.Pointer(&Am)),
        C.double(alpha), C.int(flags), C.int(N), C.int(S), C.int(L))

}

func DTrmmLower(B, A []float64, alpha float64, unit bool, ldB, ldA, N, S, L int) {
    var Bm C.mdata_t
    var Am C.mdata_t
    Bm.md =  (*C.double)(unsafe.Pointer(&B[0]))
    Bm.step = C.int(ldB)
    Am.md =  (*C.double)(unsafe.Pointer(&A[0]))
    Am.step = C.int(ldA)

    var flags Flags = LOWER
    if unit {
        flags |= UNIT
    }
    C.dmmat_trid_unb(
        (*C.mdata_t)(unsafe.Pointer(&Bm)),
        (*C.mdata_t)(unsafe.Pointer(&Am)),
        C.double(alpha), C.int(flags), C.int(N), C.int(S), C.int(L))

}

func DTrmmLowerTransA(B, A []float64, alpha float64, unit bool, ldB, ldA, N, S, L int) {
    var Bm C.mdata_t
    var Am C.mdata_t
    Bm.md =  (*C.double)(unsafe.Pointer(&B[0]))
    Bm.step = C.int(ldB)
    Am.md =  (*C.double)(unsafe.Pointer(&A[0]))
    Am.step = C.int(ldA)

    var flags Flags = LOWER|TRANSA
    if unit {
        flags |= UNIT
    }
    C.dmmat_trid_unb(
        (*C.mdata_t)(unsafe.Pointer(&Bm)),
        (*C.mdata_t)(unsafe.Pointer(&Am)),
        C.double(alpha), C.int(flags), C.int(N), C.int(S), C.int(L))

}
