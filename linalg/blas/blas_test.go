// Copyright (c) Harri Rautila, 2012

// This file is part of github.com/hrautila/linalg/blas package.
// It is free software, distributed under the terms of GNU Lesser General Public
// License Version 3, or any later version. See the COPYING tile included in this archive.

package blas

import (
	"fmt"
	"github.com/henrylee2cn/algorithm/linalg"
	"github.com/henrylee2cn/algorithm/matrix"
	"testing"
)

// a = norm2(A)
func TestDnrm2(t *testing.T) {
	fmt.Printf("* L1 * test sum: nrm2(X)\n")
	A := matrix.FloatVector([]float64{1.0, 1.0, 1.0, 1.0, 1.0, 1.0})
	v1 := Nrm2(A, &linalg.IOpt{"offset", 0})
	v2 := Nrm2(A, &linalg.IOpt{"offset", 3})
	v3 := Nrm2(A, &linalg.IOpt{"inc", 2})
	fmt.Printf("Ddnrm2\n")
	fmt.Printf("%.3f\n", v1.Float())
	fmt.Printf("%.3f\n", v2.Float())
	fmt.Printf("%.3f\n", v3.Float())
	// Output:
	// 2.499
	// 1.732
	// 1.732
}

// a = sum(X)
func TestDasum(t *testing.T) {
	fmt.Printf("* L1 * test sum: sum(X)\n")
	A := matrix.FloatVector([]float64{1.0, 1.0, 1.0, 1.0, 1.0, 1.0})
	v1 := Asum(A, &linalg.IOpt{"offset", 0})
	v2 := Asum(A, &linalg.IOpt{"offset", 3})
	v3 := Asum(A, &linalg.IOpt{"inc", 2})
	fmt.Printf("Dasum\n")
	fmt.Printf("%.3f\n", v1.Float())
	fmt.Printf("%.3f\n", v2.Float())
	fmt.Printf("%.3f\n", v3.Float())
	// Output:
	// 6.000
	// 3.000
	// 3.000
}

// v = X.T * Y
func TestDdot(t *testing.T) {
	fmt.Printf("* L1 * test dot: X.T*Y\n")
	A := matrix.FloatVector([]float64{1.0, 1.0, 1.0, 1.0, 1.0, 1.0})
	B := matrix.FloatVector([]float64{2.0, 2.0, 2.0, 2.0, 2.0, 2.0})
	v1 := Dot(A, B)
	v2 := Dot(A, B, &linalg.IOpt{"offset", 3})
	v3 := Dot(A, B, &linalg.IOpt{"inc", 2})
	fmt.Printf("Ddot: X.T * Y\n")
	fmt.Printf("%.3f\n", v1.Float())
	fmt.Printf("%.3f\n", v2.Float())
	fmt.Printf("%.3f\n", v3.Float())
	// Output:
	// 12.000
	// 6.000
	// 6.000
}

// X <--> Y
func TestDswap(t *testing.T) {
	fmt.Printf("* L1 * test swap: X <--> Y\n")
	A := matrix.FloatVector([]float64{1.0, 1.0, 1.0, 1.0, 1.0, 1.0})
	B := matrix.FloatVector([]float64{2.0, 2.0, 2.0, 2.0, 2.0, 2.0})
	Swap(A, B)
	fmt.Printf("Dswap A, B\n")
	fmt.Printf("%s\n", A)
	A = matrix.FloatVector([]float64{1.0, 1.0, 1.0, 1.0, 1.0, 1.0})
	B = matrix.FloatVector([]float64{2.0, 2.0, 2.0, 2.0, 2.0, 2.0})
	Swap(A, B, &linalg.IOpt{"offset", 3})
	fmt.Printf("Dswap A[3:], B[3:]\n")
	fmt.Printf("%s\n", A)
	A = matrix.FloatVector([]float64{1.0, 1.0, 1.0, 1.0, 1.0, 1.0})
	B = matrix.FloatVector([]float64{2.0, 2.0, 2.0, 2.0, 2.0, 2.0})
	fmt.Printf("Dswap A[::2], B[::2]\n")
	Swap(A, B, &linalg.IOpt{"inc", 2})
	fmt.Printf("%s\n", A)
}

// Dscal: X = alpha * X
func TestDscal(t *testing.T) {
	fmt.Printf("* L1 * test scal: X = alpha * X\n")
	alpha := matrix.FScalar(2.0)
	A := matrix.FloatVector([]float64{1.0, 1.0, 1.0, 1.0, 1.0, 1.0})
	Scal(A, alpha)
	fmt.Printf("Dscal 2.0 * A\n")
	fmt.Printf("%s\n", A)
	A = matrix.FloatVector([]float64{1.0, 1.0, 1.0, 1.0, 1.0, 1.0})
	Scal(A, alpha, &linalg.IOpt{"offset", 3})
	fmt.Printf("Dscal 2.0 * A[3:]\n")
	fmt.Printf("%s\n", A)
	A = matrix.FloatVector([]float64{1.0, 1.0, 1.0, 1.0, 1.0, 1.0})
	fmt.Printf("Dscal 2.0* A[::2]\n")
	Scal(A, alpha, &linalg.IOpt{"inc", 2})
	fmt.Printf("%s\n", A)
}

func TestDaxpy(t *testing.T) {
	fmt.Printf("* L1 * test axpy: Y = alpha * X + Y\n")
	X := matrix.FloatVector([]float64{1, 1, 1})
	Y := matrix.FloatVector([]float64{0, 0, 0})
	fmt.Printf("before:\nX=\n%v\nY=\n%v\n", X, Y)
	Axpy(X, Y, matrix.FScalar(5.0))
	fmt.Printf("after:\nX=\n%v\nY=\n%v\n", X, Y)
}

func TestDgemv(t *testing.T) {
	fmt.Printf("* L2 * test gemv: Y = alpha * A * X + beta * Y\n")
	A := matrix.FloatNew(3, 2, []float64{1, 1, 1, 2, 2, 2})
	X := matrix.FloatVector([]float64{1, 1})
	Y := matrix.FloatVector([]float64{0, 0, 0})
	alpha := matrix.FScalar(1.0)
	beta := matrix.FScalar(0.0)
	fmt.Printf("before: alpha=1.0, beta=0.0\nA=\n%v\nX=\n%v\nY=\n%v\n", A, X, Y)
	err := Gemv(A, X, Y, alpha, beta)
	fmt.Printf("after:\nA=\n%v\nX=\n%v\nY=\n%v\n", A, X, Y)
	fmt.Printf("* L2 * test gemv: X = alpha * A.T * Y + beta * X\n")
	err = Gemv(A, Y, X, alpha, beta, linalg.OptTrans)
	if err != nil {
		fmt.Printf("error: %s\n", err)
	}
	fmt.Printf("after:\nA=\n%v\nX=\n%v\nY=\n%v\n", A, X, Y)
}

// not run
func _TestPanic(t *testing.T) {
	PanicOnError(true)
	X := matrix.FloatWithValue(10, 1, 1.0)
	ScalFloat(X, 2.0, &linalg.IOpt{"offset", -1})
}

// Local Variables:
// tab-width: 4
// End:
