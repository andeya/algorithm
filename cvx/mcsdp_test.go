package cvx

import (
	"github.com/henrylee2cn/algorithm/cvx/sets"
	"github.com/henrylee2cn/algorithm/linalg"
	"github.com/henrylee2cn/algorithm/linalg/blas"
	"github.com/henrylee2cn/algorithm/linalg/lapack"
	"github.com/henrylee2cn/algorithm/matrix"
	"testing"
)

// Internal type for MatrixG interface
type matrixFs struct {
	n int
}

func (g *matrixFs) Gf(x, y *matrix.FloatMatrix, alpha, beta float64, trans linalg.Option) error {
	//
	// y := alpha*(-diag(x)) + beta*y
	//
	if linalg.Equal(trans, linalg.OptNoTrans) {
		blas.ScalFloat(y, beta)
		blas.AxpyFloat(x, y, -alpha, &linalg.IOpt{"incy", g.n + 1})
	} else {
		blas.ScalFloat(y, beta)
		blas.AxpyFloat(x, y, -alpha, &linalg.IOpt{"incx", g.n + 1})
	}
	return nil
}

func mcsdp(w *matrix.FloatMatrix) (*Solution, error) {
	//
	// Returns solution x, z to
	//
	//    (primal)  minimize    sum(x)
	//              subject to  w + diag(x) >= 0
	//
	//    (dual)    maximize    -tr(w*z)
	//              subject to  diag(z) = 1
	//                          z >= 0.
	//
	n := w.Rows()
	G := &matrixFs{n}

	cngrnc := func(r, x *matrix.FloatMatrix, alpha float64) (err error) {
		// Congruence transformation
		//
		//    x := alpha * r'*x*r.
		//
		// r and x are square matrices.
		//
		err = nil

		// tx = matrix(x, (n,n)) is copying and reshaping
		// scale diagonal of x by 1/2, (x is (n,n))
		tx := x.Copy()
		matrix.Reshape(tx, n, n)
		tx.Diag().Scale(0.5)

		// a := tril(x)*r
		// (python: a = +r is really making a copy of r)
		a := r.Copy()

		err = blas.TrmmFloat(tx, a, 1.0, linalg.OptLeft)

		// x := alpha*(a*r' + r*a')
		err = blas.Syr2kFloat(r, a, tx, alpha, 0.0, linalg.OptTrans)

		// x[:] = tx[:]
		tx.CopyTo(x)
		return
	}

	Fkkt := func(W *sets.FloatMatrixSet) (KKTFunc, error) {

		//    Solve
		//                  -diag(z)                           = bx
		//        -diag(x) - inv(rti*rti') * z * inv(rti*rti') = bs
		//
		//    On entry, x and z contain bx and bs.
		//    On exit, they contain the solution, with z scaled
		//    (inv(rti)'*z*inv(rti) is returned instead of z).
		//
		//    We first solve
		//
		//        ((rti*rti') .* (rti*rti')) * x = bx - diag(t*bs*t)
		//
		//    and take z  = -rti' * (diag(x) + bs) * rti.

		var err error = nil
		rti := W.At("rti")[0]

		// t = rti*rti' as a nonsymmetric matrix.
		t := matrix.FloatZeros(n, n)
		err = blas.GemmFloat(rti, rti, t, 1.0, 0.0, linalg.OptTransB)
		if err != nil {
			return nil, err
		}

		// Cholesky factorization of tsq = t.*t.
		tsq := matrix.Mul(t, t)
		err = lapack.Potrf(tsq)
		if err != nil {
			return nil, err
		}

		f := func(x, y, z *matrix.FloatMatrix) (err error) {
			// tbst := t * zs * t = t * bs * t
			tbst := z.Copy()
			matrix.Reshape(tbst, n, n)
			cngrnc(t, tbst, 1.0)

			// x := x - diag(tbst) = bx - diag(rti*rti' * bs * rti*rti')
			diag := tbst.Diag().Transpose()
			x.Minus(diag)

			// x := (t.*t)^{-1} * x = (t.*t)^{-1} * (bx - diag(t*bs*t))
			err = lapack.Potrs(tsq, x)

			// z := z + diag(x) = bs + diag(x)
			// z, x are really column vectors here
			z.AddIndexes(matrix.MakeIndexSet(0, n*n, n+1), x.FloatArray())

			// z := -rti' * z * rti = -rti' * (diag(x) + bs) * rti
			cngrnc(rti, z, -1.0)
			return nil
		}
		return f, nil
	}

	c := matrix.FloatWithValue(n, 1, 1.0)

	// initial feasible x: x = 1.0 - min(lmbda(w))
	lmbda := matrix.FloatZeros(n, 1)
	wp := w.Copy()
	lapack.Syevx(wp, lmbda, nil, 0.0, nil, []int{1, 1}, linalg.OptRangeInt)
	x0 := matrix.FloatZeros(n, 1).Add(-lmbda.GetAt(0, 0) + 1.0)
	s0 := w.Copy()
	s0.Diag().Plus(x0.Transpose())
	matrix.Reshape(s0, n*n, 1)

	// initial feasible z is identity
	z0 := matrix.FloatIdentity(n)
	matrix.Reshape(z0, n*n, 1)

	dims := sets.DSetNew("l", "q", "s")
	dims.Set("s", []int{n})

	primalstart := sets.FloatSetNew("x", "s")
	dualstart := sets.FloatSetNew("z")
	primalstart.Set("x", x0)
	primalstart.Set("s", s0)
	dualstart.Set("z", z0)

	var solopts SolverOptions
	solopts.MaxIter = 30
	solopts.ShowProgress = false
	h := w.Copy()
	matrix.Reshape(h, h.NumElements(), 1)
	return ConeLpCustomMatrix(c, G, h, nil, nil, dims, Fkkt, &solopts, primalstart, dualstart)
}

func TestMcSdp(t *testing.T) {
	dataref := []float64{
		0.13391860811867590, -0.08810099183143839, 1.67440840625377385,
		0.73364110729257948, 0.99752463160201232, -1.27750208100276641,
		-2.39671528273320744, -0.67928016472921104, -0.03909131843358723,
		0.89355554552081806, -0.01764779458978329, -1.29655690877732033,
		-0.66798050600863512, 0.18170877830799323, 0.83105079034069296,
		-0.54824847139682842, -0.63803150269339892, 0.00708853169542829,
		-0.66985369437758802, -0.82776233138694755, 0.61250003036052025,
		-0.36758841015776106, -0.28949340911906524, -0.91687467081464402,
		-0.58380545879822354, 0.06461529015869494, 0.04674385977341894,
		-0.75072697866800864, -0.14423268736469502, 1.47556079692186803,
		-0.20647897491996936, -0.26069727887528138, -0.13821797189084586,
		-2.54339721705159461, 0.57741853670968823, -0.03610986310676272,
		-0.22617256815925429, 0.71500221712914425, 1.57726688654131597,
		-0.39227386148741172, 1.18277085579223740, 1.79311368741147459,
		-0.73134596730591273, -0.62981768127016369, 2.07069457080308439,
		-0.90300600743388137, 1.73554553761132668, 1.38308886012672638,
		-0.37536663987692281, -0.30431019951136606, 1.25357124882638882,
		-0.27384517924722629, -0.04632938142394107, 0.79898247010790735,
		0.53995588924910287, 0.08585033894586751, 1.35529122074053832,
		0.01385086083133875, 0.75136369502320788, -1.57490917820568810,
		-0.23737825305680132, 0.33649059262790887, -0.07959967279201274,
		-0.54087224229887154, -1.38043394277880926, 1.71655186027439033,
		0.41811180400013570, 1.52495658744607199, 0.68647633127251673,
		-0.58855597085742983, 1.31674790938058739, -0.83298414482940919,
		0.77396029737830319, -0.80351540646387032, -0.08543027502251639,
		1.49619789524058855, 0.39776679351920241, 3.33842986686760090,
		1.80093502635427316, 1.75627861060973078, 0.66747365962131211,
		-1.05016810484723888, 0.32937053624566975, 0.45508103845913089,
		-2.34700732048193306, -0.61389923742031549, -2.16180058868422442,
		-1.00134312023216387, -0.46226762826081341, -0.81756106137679618,
		0.18439203820434547, -0.39491697602688680, -1.66955730897081911,
		-0.94514784362504667, 0.00356999348571242, 1.15573946640829073,
		1.87743135406613315, 0.11219893028802071, 0.18121323774935597,
		1.27830445029226136}

	xref := []float64{
		4.35021720027799841, 3.05093154317696458, 3.81967245723734372,
		4.46546919954538968, 1.84870016923725555, 1.98581804497831516,
		4.10446362983979540, 0.71734387468946914, 4.36803763631339592, 4.18804903479807233}

	zref := []float64{
		1.00000000000000022, 0.85303260971097450, -0.85020234203223444,
		-0.64630863726711674, -0.66211228030949409, 0.12515757728046226,
		0.96725545921055622, -0.37698782755909727, 0.39233601196202428,
		0.43202712322170089, 0.85303260971097450, 0.99999999999999956,
		-0.99998464240308171, -0.94953899113609297, -0.17372147858022127,
		0.62451665675716073, 0.69265128272781040, -0.80493649049792482,
		0.81469110300198433, 0.83917543899780545, -0.85020234203223444,
		-0.99998464240308171, 1.00000000000000089, 0.95121850038903877,
		0.16840142730281188, -0.62872513058587987, -0.68874630140860726,
		0.80812857148627615, -0.81781005005790097, -0.84210005596350634,
		-0.64630863726711674, -0.94953899113609297, 0.95121850038903877,
		1.00000000000000000, -0.14392344681149352, -0.83796540396053343,
		-0.43147384025661650, 0.95042502132018603, -0.95546406276523821,
		-0.96741055070299309, -0.66211228030949409, -0.17372147858022127,
		0.16840142730281188, -0.14392344681149352, 1.00000000000000067,
		0.66064328265023031, -0.83063373455540346, -0.44450355615883147,
		0.42954780821049410, 0.38980770530806869, 0.12515757728046226,
		0.62451665675716073, -0.62872513058587987, -0.83796540396053343,
		0.66064328265023031, 0.99999999999999978, -0.13074907207216546,
		-0.96611746250486030, 0.96169189357427165, 0.94884018282016180,
		0.96725545921055622, 0.69265128272781040, -0.68874630140860726,
		-0.43147384025661650, -0.83063373455540346, -0.13074907207216546,
		1.00000000000000044, -0.12956585354159514, 0.14603486320665809,
		0.18898493355529175, -0.37698782755909727, -0.80493649049792482,
		0.80812857148627615, 0.95042502132018603, -0.44450355615883147,
		-0.96611746250486030, -0.12956585354159514, 1.00000000000000067,
		-0.99986146360473316, -0.99818853872682889, 0.39233601196202428,
		0.81469110300198433, -0.81781005005790097, -0.95546406276523821,
		0.42954780821049410, 0.96169189357427165, 0.14603486320665809,
		-0.99986146360473316, 1.00000000000000022, 0.99905052020024310,
		0.43202712322170089, 0.83917543899780545, -0.84210005596350634,
		-0.96741055070299309, 0.38980770530806869, 0.94884018282016180,
		0.18898493355529175, -0.99818853872682889, 0.99905052020024310, 0.99999999999999933}

	data := matrix.FloatNew(10, 10, dataref)

	sol, err := mcsdp(data)
	if sol != nil && sol.Status == Optimal {
		fail := false
		x := sol.Result.At("x")[0]
		z := sol.Result.At("z")[0]
		//matrix.Reshape(z, data.Rows(), data.Rows())
		t.Logf("x=\n%v\n", x.ToString("%.9f"))
		xe, _ := nrmError(matrix.FloatVector(xref), x)
		if xe > TOL {
			t.Logf("x differs [%.3e] from exepted too much.", xe)
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
		t.Logf("status: %v\n", err)
		t.Fail()
	}
}

// Local Variables:
// tab-width: 4
// End:
