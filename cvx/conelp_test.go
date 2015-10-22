package cvx

import (
	"github.com/henrylee2cn/algorithm/cvx/sets"
	"github.com/henrylee2cn/algorithm/linalg/blas"
	"github.com/henrylee2cn/algorithm/matrix"
	"testing"
)

const TOL = 1e-8

func nrmError(ref, val *matrix.FloatMatrix) (nrm float64, diff *matrix.FloatMatrix) {
	diff = ref.Minus(val)
	nrm = blas.Nrm2(diff).Float()
	return
}

func TestConeLp(t *testing.T) {

	gdata := [][]float64{
		[]float64{16., 7., 24., -8., 8., -1., 0., -1., 0., 0., 7.,
			-5., 1., -5., 1., -7., 1., -7., -4.},
		[]float64{-14., 2., 7., -13., -18., 3., 0., 0., -1., 0., 3.,
			13., -6., 13., 12., -10., -6., -10., -28.},
		[]float64{5., 0., -15., 12., -6., 17., 0., 0., 0., -1., 9.,
			6., -6., 6., -7., -7., -6., -7., -11.}}

	hdata := []float64{-3., 5., 12., -2., -14., -13., 10., 0., 0., 0., 68.,
		-30., -19., -30., 99., 23., -19., 23., 10.}

	// these reference values obtained from running cvxopt conelp.py example
	xref := []float64{-1.22091525026262993, 0.09663323966626469, 3.57750155386611057}

	sref := []float64{
		0.00000172588537019, 13.35314040819201864,
		94.28805677232460880, -53.44110853283719109,
		18.97172963929198275, -75.32834138499130461,
		10.00000013568614321, -1.22091525026262993,
		0.09663323966626476, 3.57750155386611146,
		44.05899318373081286, -58.82581769017131990,
		4.26572401145687596, -58.82581769017131990,
		124.10382738701650851, 40.46243652188705653,
		4.26572401145687596, 40.46243652188705653,
		47.17458693781828316}

	zref := []float64{
		0.09299833991484617, 0.00000001060210894,
		0.23532251654806322, 0.13337937743566930,
		-0.04734875722474355, 0.18800192060450249,
		0.00000001245876667, 0.00000000007816348,
		-0.00000000039584268, -0.00000000183463577,
		0.12558704894101563, 0.08777794737598217,
		-0.08664401207348003, 0.08777794737598217,
		0.06135161787371416, -0.06055906182304811,
		-0.08664401207348003, -0.06055906182304811,
		0.05977675078191153}

	c := matrix.FloatVector([]float64{-6., -4., -5.})
	G := matrix.FloatMatrixFromTable(gdata)
	h := matrix.FloatVector(hdata)

	dims := sets.DSetNew("l", "q", "s")
	dims.Set("l", []int{2})
	dims.Set("q", []int{4, 4})
	dims.Set("s", []int{3})

	var solopts SolverOptions
	solopts.MaxIter = 30
	solopts.ShowProgress = false
	sol, err := ConeLp(c, G, h, nil, nil, dims, &solopts, nil, nil)
	if err == nil {
		fail := false
		x := sol.Result.At("x")[0]
		s := sol.Result.At("s")[0]
		z := sol.Result.At("z")[0]
		t.Logf("Optimal\n")
		t.Logf("x=\n%v\n", x.ToString("%.9f"))
		t.Logf("s=\n%v\n", s.ToString("%.9f"))
		t.Logf("z=\n%v\n", z.ToString("%.9f"))
		xe, _ := nrmError(matrix.FloatVector(xref), x)
		if xe > TOL {
			t.Logf("x differs [%.3e] from exepted too much.", xe)
			fail = true
		}
		se, _ := nrmError(matrix.FloatVector(sref), s)
		if se > TOL {
			t.Logf("s differs [%.3e] from exepted too much.", se)
			fail = true
		}
		ze, _ := nrmError(matrix.FloatVector(zref), z)
		if ze > TOL {
			t.Logf("z differs [%.3e] from exepted too much.", ze)
			fail = true
		}
		if fail {
			t.Fail()
		}
	} else {
		t.Logf("status: %s\n", err)
		t.Fail()
	}
}

// Local Variables:
// tab-width: 4
// End:
