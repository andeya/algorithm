// Copyright (c) Harri Rautila, 2012

// This file is part of  package. It is free software, distributed
// under the terms of GNU Lesser General Public License Version 3, or any later
// version. See the COPYING tile included in this archive.

package lapack

// Analytic centering example at the end of chapter 4 CVXOPT package documentation.
// This file adaptation  of the original acent.py example.

import (
	"errors"
	"fmt"
	"github.com/henrylee2cn/algorithm/linalg"
	"github.com/henrylee2cn/algorithm/linalg/blas"
	"github.com/henrylee2cn/algorithm/matrix"
	"math"
	"testing"
)

const (
	MAXITERS = 100
	ALPHA    = 0.01
	BETA     = 0.5
	TOL      = 1e-8
)

// Computes analytic center of A*x <= b with A m by n of rank n.
// We assume that b > 0 and the feasible set is bounded.
func acent(A, b *matrix.FloatMatrix, niters int) (x *matrix.FloatMatrix, ntdecrs []float64, err error) {

	err = nil
	if niters <= 0 {
		niters = MAXITERS
	}
	ntdecrs = make([]float64, 0, niters)

	if A.Rows() != b.Rows() {
		return nil, nil, errors.New("A.Rows() != b.Rows()")
	}

	m, n := A.Size()
	x = matrix.FloatZeros(n, 1)
	H := matrix.FloatZeros(n, n)
	// Helper m*n matrix
	Dmn := matrix.FloatZeros(m, n)

	for i := 0; i < niters; i++ {

		// Gradient is g = A^T * (1.0/(b - A*x)). d = 1.0/(b - A*x)
		// d is m*1 matrix, g is n*1 matrix
		d := matrix.Minus(b, matrix.Times(A, x)).Inv()
		g := matrix.Times(A.Transpose(), d)

		// Hessian is H = A^T * diag(1./(b-A*x))^2 * A.
		// in the original python code expression d[:,n*[0]] creates
		// a m*n matrix where each column is copy of column 0.
		// We do it here manually.
		for i := 0; i < n; i++ {
			Dmn.SetColumn(i, d)
		}

		// Function mul creates element wise product of matrices.
		Asc := matrix.Mul(Dmn, A)
		blas.SyrkFloat(Asc, H, 1.0, 0.0, linalg.OptTrans)

		// Newton step is v = H^-1 * g.
		v := matrix.Scale(g, -1.0)
		PosvFloat(H, v)

		// Directional derivative and Newton decrement.
		lam := blas.DotFloat(g, v)
		ntdecrs = append(ntdecrs, math.Sqrt(-lam))
		if ntdecrs[len(ntdecrs)-1] < TOL {
			return x, ntdecrs, err
		}

		// Backtracking line search.
		// y = d .* A*v
		y := matrix.Mul(d, matrix.Times(A, v))
		step := 1.0
		for 1-step*y.Max() < 0 {
			step *= BETA
		}

	search:
		for {
			// t = -step*y + 1 [e.g. t = 1 - step*y]
			t := matrix.Scale(y, -step).Add(1.0)

			// ts = sum(log(1-step*y))
			ts := t.Log().Sum()
			if -ts < ALPHA*step*lam {
				break search
			}
			step *= BETA
		}
		v.Scale(step)
		x.Plus(v)
	}
	// no solution !!
	err = errors.New(fmt.Sprintf("Iteration %d exhausted\n", niters))
	return x, ntdecrs, err
}

func TestAcent(t *testing.T) {
	// matrix string in row order presentation
	Adata := [][]float64{
		[]float64{-7.44e-01, 1.11e-01, 1.29e+00, 2.62e+00, -1.82e+00},
		[]float64{4.59e-01, 7.06e-01, 3.16e-01, -1.06e-01, 7.80e-01},
		[]float64{-2.95e-02, -2.22e-01, -2.07e-01, -9.11e-01, -3.92e-01},
		[]float64{-7.75e-01, 1.03e-01, -1.22e+00, -5.74e-01, -3.32e-01},
		[]float64{-1.80e+00, 1.24e+00, -2.61e+00, -9.31e-01, -6.38e-01}}

	bdata := []float64{
		8.38e-01, 9.92e-01, 9.56e-01, 6.14e-01, 6.56e-01,
		3.57e-01, 6.36e-01, 5.08e-01, 8.81e-03, 7.08e-02}

	// these are solution obtained from running cvxopt acent.py with above data
	solData := []float64{-11.59728373909344512, -1.35196389161339936,
		7.21894899350256303, -3.29159917142051528, 4.90454147385329176}

	ntData := []float64{
		1.5163484265903457, 1.2433928210771914, 1.0562922103520955, 0.8816246051011607,
		0.7271128861543598, 0.42725003346248974, 0.0816777301914883, 0.0005458037072843131,
		1.6259980735305693e-10}

	b := matrix.FloatVector(bdata)
	Al := matrix.FloatMatrixFromTable(Adata, matrix.RowOrder)
	Au := matrix.Scale(Al, -1.0)
	A := matrix.FloatZeros(2*Al.Rows(), Al.Cols())
	A.SetSubMatrix(0, 0, Al)
	A.SetSubMatrix(Al.Rows(), 0, Au)

	x, nt, err := acent(A, b, 10)
	if err != nil {
		t.Logf("Acent error: %s", err)
		t.Fail()
	}
	solref := matrix.FloatVector(solData)
	ntref := matrix.FloatVector(ntData)
	soldf := matrix.Minus(x, solref)
	ntdf := matrix.Minus(matrix.FloatVector(nt), ntref)
	solNrm := blas.Nrm2Float(soldf)
	ntNrm := blas.Nrm2Float(ntdf)
	t.Logf("x  [diff=%.2e]:\n%v\n", solNrm, x)
	t.Logf("nt [diff=%.2e]:\n%v\n", ntNrm, nt)

	if solNrm > TOL {
		t.Log("solution deviates too much from expected\n")
		t.Fail()
	}
}

// Local Variables:
// tab-width: 4
// End:
