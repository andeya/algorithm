// This file is part of cvx package. It is free software, distributed
// under the terms of GNU Lesser General Public License Version 3, or any later
// version. See the COPYING tile included in this archive.

package cvx

import (
	"github.com/henrylee2cn/algorithm/cvx/sets"
	"github.com/henrylee2cn/algorithm/linalg/blas"
	"github.com/henrylee2cn/algorithm/matrix"
	"testing"
)

func derror(ref, val *matrix.FloatMatrix) (nrm float64, diff *matrix.FloatMatrix) {
	diff = ref.Minus(val)
	nrm = blas.Nrm2(diff).Float()
	return
}

func solve(solver string, c, G, h, A, b *matrix.FloatMatrix) (*Solution, error) {

	var solopts SolverOptions
	solopts.MaxIter = 30
	solopts.ShowProgress = true
	solopts.KKTSolverName = solver
	return Lp(c, G, h, A, b, &solopts, nil, nil)
}

// This simple test will succeed on initial points, no iterations.

// min(x) st. x+y = 1, x >= y
func TestSimple(t *testing.T) {

	A := matrix.FloatNew(2, 3, []float64{1.0, -1.0, 0.0, 1.0, 0.0, 1.0})
	b := matrix.FloatNew(2, 1, []float64{1.0, 0.0})
	c := matrix.FloatNew(3, 1, []float64{0.0, 1.0, 0.0})
	G := matrix.FloatNew(1, 3, []float64{0.0, -1.0, 1.0})
	h := matrix.FloatNew(1, 1, []float64{0.0})
	dims := sets.NewDimensionSet("l", "q", "s")
	dims.Set("l", []int{1})

	t.Logf("A=\n%v\n", A)
	t.Logf("b=\n%v\n", b)
	t.Logf("G=\n%v\n", G)
	t.Logf("h=\n%v\n", h)
	t.Logf("c=\n%v\n", c)

	// this should work...
	t.Logf("Ldl solver ...\n")
	sol, err := solve("ldl", c, G, h, A, b)
	if sol != nil && sol.Status == Optimal {
		x := sol.Result.At("x")[0]
		s := sol.Result.At("s")[0]
		z := sol.Result.At("z")[0]
		t.Logf("x=\n%v\n", x.ToString("%.9f"))
		t.Logf("s=\n%v\n", s.ToString("%.9f"))
		t.Logf("z=\n%v\n", z.ToString("%.9f"))
	} else {
		t.Logf("status: %v\n", err)
		t.Fail()
	}

	// this should work too
	t.Logf("chol2 solver ...\n")
	sol, err = solve("chol2", c, G, h, A, b)
	if err != nil {
		t.Logf("chol2 status: %v\n", err)
		t.Fail()
	}
}

// Local Variables:
// tab-width: 4
// indent-tabs-mode: nil
// End:
