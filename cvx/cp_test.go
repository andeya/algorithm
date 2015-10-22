// Copyright (c) Harri Rautila, 2012

// This file is part of github.com/hrautila/cvx package.
// It is free software, distributed under the terms of GNU Lesser General Public
// License Version 3, or any later version. See the COPYING tile included in this archive.

package cvx

import (
	"github.com/henrylee2cn/algorithm/matrix"
	//"github.com/henrylee2cn/algorithm/linalg/blas"
	"errors"
	"github.com/henrylee2cn/algorithm/cvx/sets"
	"testing"
)

type acenterProg struct {
	rows, cols int
}

func (p *acenterProg) F0() (mnl int, x0 *matrix.FloatMatrix, err error) {
	err = nil
	mnl = 0
	x0 = matrix.FloatZeros(p.rows, p.cols)
	return
}

func (p *acenterProg) F1(x *matrix.FloatMatrix) (f, Df *matrix.FloatMatrix, err error) {
	f = nil
	Df = nil
	err = nil
	max := matrix.Abs(x).Max()
	if max >= 1.0 {
		err = errors.New("max(abs(x)) >= 1.0")
		return
	}
	// u = 1 - x**2
	u := matrix.Pow(x, 2.0).Scale(-1.0).Add(1.0)
	val := -matrix.Log(u).Sum()
	f = matrix.FloatValue(val)
	Df = matrix.Div(matrix.Scale(x, 2.0), u).Transpose()
	return
}

func (p *acenterProg) F2(x, z *matrix.FloatMatrix) (f, Df, H *matrix.FloatMatrix, err error) {
	f, Df, err = p.F1(x)
	u := matrix.Pow(x, 2.0).Scale(-1.0).Add(1.0)
	z0 := z.GetIndex(0)
	u2 := matrix.Pow(u, 2.0)
	hd := matrix.Div(matrix.Add(u2, 1.0), u2).Scale(2 * z0)
	H = matrix.FloatDiagonal(hd.NumElements(), hd.FloatArray()...)
	return
}

// The analytic centering with cone constraints example of section 9.1
// (Problems with nonlinear objectives).
func TestCp(t *testing.T) {

	xref := []float64{0.41132359189354400, 0.55884774432611484, -0.72007090016957931}

	F := &acenterProg{3, 1}

	gdata := [][]float64{
		[]float64{0., -1., 0., 0., -21., -11., 0., -11., 10., 8., 0., 8., 5.},
		[]float64{0., 0., -1., 0., 0., 10., 16., 10., -10., -10., 16., -10., 3.},
		[]float64{0., 0., 0., -1., -5., 2., -17., 2., -6., 8., -17., -7., 6.}}

	G := matrix.FloatMatrixFromTable(gdata)
	h := matrix.FloatVector(
		[]float64{1.0, 0.0, 0.0, 0.0, 20., 10., 40., 10., 80., 10., 40., 10., 15.})

	var solopts SolverOptions
	solopts.MaxIter = 40
	solopts.ShowProgress = false

	dims := sets.NewDimensionSet("l", "q", "s")
	dims.Set("l", []int{0})
	dims.Set("q", []int{4})
	dims.Set("s", []int{3})

	sol, err := Cp(F, G, h, nil, nil, dims, &solopts)
	if err == nil && sol.Status == Optimal {
		x := sol.Result.At("x")[0]
		t.Logf("x = \n%v\n", x.ToString("%.9f"))
		xe, _ := nrmError(matrix.FloatVector(xref), x)
		if xe > TOL {
			t.Logf("x differs [%.3e] from exepted too much.", xe)
			t.Fail()
		}
	} else {
		t.Logf("result: %v\n", err)
		t.Fail()
	}
}

// Local Variables:
// tab-width: 4
// End:
