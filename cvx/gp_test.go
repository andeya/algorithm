package cvx

import (
	"github.com/henrylee2cn/algorithm/matrix"
	"testing"
)

// The small GP of section 9.3 (Geometric programming).
func TestGp(t *testing.T) {

	xref := []float64{1.06032641296944741, 1.75347359157296845, 2.44603683900611868}

	aflr := 1000.0
	awall := 100.0
	alpha := 0.5
	beta := 2.0
	gamma := 0.5
	delta := 2.0

	fdata := [][]float64{
		[]float64{-1.0, 1.0, 1.0, 0.0, -1.0, 1.0, 0.0, 0.0},
		[]float64{-1.0, 1.0, 0.0, 1.0, 1.0, -1.0, 1.0, -1.0},
		[]float64{-1.0, 0.0, 1.0, 1.0, 0.0, 0.0, -1.0, 1.0}}

	gdata := []float64{1.0, 2.0 / awall, 2.0 / awall, 1.0 / aflr, alpha, 1.0 / beta, gamma, 1.0 / delta}

	g := matrix.FloatNew(8, 1, gdata).Log()
	F := matrix.FloatMatrixFromTable(fdata)
	K := []int{1, 2, 1, 1, 1, 1, 1}

	var solopts SolverOptions
	solopts.MaxIter = 40
	solopts.ShowProgress = false
	solopts.KKTSolverName = "ldl"
	sol, err := Gp(K, F, g, nil, nil, nil, nil, &solopts)
	if sol != nil && sol.Status == Optimal {
		x := sol.Result.At("x")[0]
		r := matrix.Exp(x)
		h := r.GetIndex(0)
		w := r.GetIndex(1)
		d := r.GetIndex(2)
		t.Logf("x=\n%v\n", x.ToString("%.9f"))
		t.Logf("h = %f,  w = %f, d = %f.\n", h, w, d)
		xe, _ := nrmError(matrix.FloatVector(xref), x)
		if xe > TOL {
			t.Logf("x differs [%.3e] from exepted too much.", xe)
			t.Fail()
		}
	} else {
		t.Logf("status: %v\n", err)
		t.Fail()
	}
}

// Local Variables:
// tab-width: 4
// End:
