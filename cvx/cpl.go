// Copyright (c) Harri Rautila, 2012

// This file is part of github.com/hrautila/cvx package.
// It is free software, distributed under the terms of GNU Lesser General Public
// License Version 3, or any later version. See the COPYING tile included in this archive.

package cvx

import (
	"errors"
	"fmt"
	"github.com/henrylee2cn/algorithm/cvx/checkpnt"
	"github.com/henrylee2cn/algorithm/cvx/sets"
	la "github.com/henrylee2cn/algorithm/linalg"
	"github.com/henrylee2cn/algorithm/linalg/blas"
	"github.com/henrylee2cn/algorithm/matrix"
	"math"
)

// Solves a convex optimization problem with a linear objective
//
//        minimize    c'*x
//        subject to  f(x) <= 0
//                    G*x <= h
//                    A*x = b.
//
// f is vector valued, convex and twice differentiable.  The linear
// inequalities are with respect to a cone C defined as the Cartesian
// product of N + M + 1 cones:
//
//        C = C_0 x C_1 x .... x C_N x C_{N+1} x ... x C_{N+M}.
//
// The first cone C_0 is the nonnegative orthant of dimension ml.  The
// next N cones are second order cones of dimension r[0], ..., r[N-1].
// The second order cone of dimension m is defined as
//
//        { (u0, u1) in R x R^{m-1} | u0 >= ||u1||_2 }.
//
// The next M cones are positive semidefinite cones of order t[0], ..., t[M-1] >= 0.
//
// The structure of C is specified by DimensionSet dims which holds following sets
//
//   dims.At("l")  l, the dimension of the nonnegative orthant (array of length 1)
//   dims.At("q")  r[0], ... r[N-1], list with the dimesions of the second-order cones
//   dims.At("s")  t[0], ... t[M-1], array with the dimensions of the positive
//                 semidefinite cones
//
// The default value for dims is l: []int{h.Rows()}, q: []int{}, s: []int{}.
//
// On exit Solution contains the result and information about the accurancy of the
// solution. if SolutionStatus is Optimal then Solution.Result contains solutions
// for the problems.
//
//   Result.At("x")[0]    primal solution
//   Result.At("snl")[0]  non-linear constraint slacks
//   Result.At("sl")[0]   linear constraint slacks
//   Result.At("y")[0]    values for linear equality constraints y
//   Result.At("znl")[0]  values of dual variables for nonlinear inequalities
//   Result.At("zl")[0]   values of dual variables for linear inequalities
//
// If err is non-nil then sol is nil and err contains information about the argument or
// computation error.
//
func Cpl(F ConvexProg, c, G, h, A, b *matrix.FloatMatrix, dims *sets.DimensionSet, solopts *SolverOptions) (sol *Solution, err error) {

	var mnl int
	var x0 *matrix.FloatMatrix

	mnl, x0, err = F.F0()
	if err != nil {
		return
	}

	if x0.Cols() != 1 {
		err = errors.New("'x0' must be matrix with one column")
		return
	}
	if c == nil {
		err = errors.New("'c' must be non nil matrix")
		return
	}
	if !c.SizeMatch(x0.Size()) {
		err = errors.New(fmt.Sprintf("'c' must be matrix of size (%d,1)", x0.Rows()))
		return
	}

	if h == nil {
		h = matrix.FloatZeros(0, 1)
	}
	if h.Cols() > 1 {
		err = errors.New("'h' must be matrix with 1 column")
		return
	}

	if dims == nil {
		dims = sets.NewDimensionSet("l", "q", "s")
		dims.Set("l", []int{h.Rows()})
	}

	cdim := dims.Sum("l", "q") + dims.SumSquared("s")
	//cdim_pckd := dims.Sum("l", "q") + dims.SumPacked("s")
	//cdim_diag := dims.Sum("l", "q", "s")

	if h.Rows() != cdim {
		err = errors.New(fmt.Sprintf("'h' must be float matrix of size (%d,1)", cdim))
		return
	}

	if G == nil {
		G = matrix.FloatZeros(0, c.Rows())
	}
	if !G.SizeMatch(cdim, c.Rows()) {
		estr := fmt.Sprintf("'G' must be of size (%d,%d)", cdim, c.Rows())
		err = errors.New(estr)
		return
	}

	// Check A and set defaults if it is nil
	if A == nil {
		// zeros rows reduces Gemv to vector products
		A = matrix.FloatZeros(0, c.Rows())
	}
	if A.Cols() != c.Rows() {
		estr := fmt.Sprintf("'A' must have %d columns", c.Rows())
		err = errors.New(estr)
		return
	}

	// Check b and set defaults if it is nil
	if b == nil {
		b = matrix.FloatZeros(0, 1)
	}
	if b.Cols() != 1 {
		estr := fmt.Sprintf("'b' must be a matrix with 1 column")
		err = errors.New(estr)
		return
	}
	if b.Rows() != A.Rows() {
		estr := fmt.Sprintf("'b' must have length %d", A.Rows())
		err = errors.New(estr)
		return
	}

	var mc = matrixVar{c}
	var mb = matrixVar{b}
	var mA = matrixVarA{A}
	var mG = matrixVarG{G, dims}

	solvername := solopts.KKTSolverName
	if len(solvername) == 0 {
		if len(dims.At("q")) > 0 || len(dims.At("s")) > 0 {
			solvername = "chol"
		} else {
			solvername = "chol2"
		}
	}

	var factor kktFactor
	var kktsolver KKTCpSolver = nil
	if kktfunc, ok := solvers[solvername]; ok {
		// kkt function returns us problem spesific factor function.
		factor, err = kktfunc(G, dims, A, mnl)
		// solver is
		kktsolver = func(W *sets.FloatMatrixSet, x, z *matrix.FloatMatrix) (KKTFunc, error) {
			_, Df, H, err := F.F2(x, z)
			if err != nil {
				return nil, err
			}
			return factor(W, H, Df)
		}
	} else {
		err = errors.New(fmt.Sprintf("solver '%s' not known", solvername))
		return
	}

	//return CplCustom(F, c, &mG, h, &mA, b, dims, kktsolver, solopts)
	return cpl_problem(F, &mc, &mG, h, &mA, &mb, dims, kktsolver, solopts, x0, mnl)
}

// Solves a convex optimization problem with a linear objective
//
//        minimize    c'*x
//        subject to  f(x) <= 0
//                    G*x <= h
//                    A*x = b.
//
// using custom KTT equation solver.
//
func CplCustomKKT(F ConvexProg, c *matrix.FloatMatrix, G, h, A, b *matrix.FloatMatrix,
	dims *sets.DimensionSet, kktsolver KKTCpSolver,
	solopts *SolverOptions) (sol *Solution, err error) {

	var mnl int
	var x0 *matrix.FloatMatrix

	mnl, x0, err = F.F0()
	if err != nil {
		return
	}

	if x0.Cols() != 1 {
		err = errors.New("'x0' must be matrix with one column")
		return
	}
	if c == nil {
		err = errors.New("'c' must be non nil matrix")
		return
	}
	if !c.SizeMatch(x0.Size()) {
		err = errors.New(fmt.Sprintf("'c' must be matrix of size (%d,1)", x0.Rows()))
		return
	}

	if h == nil {
		h = matrix.FloatZeros(0, 1)
	}
	if h.Cols() > 1 {
		err = errors.New("'h' must be matrix with 1 column")
		return
	}

	if dims == nil {
		dims = sets.NewDimensionSet("l", "q", "s")
		dims.Set("l", []int{h.Rows()})
	}

	cdim := dims.Sum("l", "q") + dims.SumSquared("s")

	if h.Rows() != cdim {
		err = errors.New(fmt.Sprintf("'h' must be float matrix of size (%d,1)", cdim))
		return
	}

	if G == nil {
		G = matrix.FloatZeros(0, c.Rows())
	}
	if !G.SizeMatch(cdim, c.Rows()) {
		estr := fmt.Sprintf("'G' must be of size (%d,%d)", cdim, c.Rows())
		err = errors.New(estr)
		return
	}

	// Check A and set defaults if it is nil
	if A == nil {
		// zeros rows reduces Gemv to vector products
		A = matrix.FloatZeros(0, c.Rows())
	}
	if A.Cols() != c.Rows() {
		estr := fmt.Sprintf("'A' must have %d columns", c.Rows())
		err = errors.New(estr)
		return
	}

	// Check b and set defaults if it is nil
	if b == nil {
		b = matrix.FloatZeros(0, 1)
	}
	if b.Cols() != 1 {
		estr := fmt.Sprintf("'b' must be a matrix with 1 column")
		err = errors.New(estr)
		return
	}
	if b.Rows() != A.Rows() {
		estr := fmt.Sprintf("'b' must have length %d", A.Rows())
		err = errors.New(estr)
		return
	}

	var mc = matrixVar{c}
	var mb = matrixVar{b}
	var mA = matrixVarA{A}
	var mG = matrixVarG{G, dims}

	return cpl_problem(F, &mc, &mG, h, &mA, &mb, dims, kktsolver, solopts, x0, mnl)
}

// Solves a convex optimization problem with a linear objective
//
//        minimize    c'*x
//        subject to  f(x) <= 0
//                    G*x <= h
//                    A*x = b.
//
// using custom KTT equation solver and custom constraints G and A.
//
func CplCustomMatrix(F ConvexProg, c *matrix.FloatMatrix, G MatrixG, h *matrix.FloatMatrix,
	A MatrixA, b *matrix.FloatMatrix, dims *sets.DimensionSet, kktsolver KKTCpSolver,
	solopts *SolverOptions) (sol *Solution, err error) {

	var mnl int
	var x0 *matrix.FloatMatrix

	mnl, x0, err = F.F0()
	if err != nil {
		return
	}

	if x0.Cols() != 1 {
		err = errors.New("'x0' must be matrix with one column")
		return
	}
	if c == nil {
		err = errors.New("'c' must be non nil matrix")
		return
	}
	if !c.SizeMatch(x0.Size()) {
		err = errors.New(fmt.Sprintf("'c' must be matrix of size (%d,1)", x0.Rows()))
		return
	}

	if h == nil {
		h = matrix.FloatZeros(0, 1)
	}
	if h.Cols() > 1 {
		err = errors.New("'h' must be matrix with 1 column")
		return
	}

	if dims == nil {
		dims = sets.NewDimensionSet("l", "q", "s")
		dims.Set("l", []int{h.Rows()})
	}

	cdim := dims.Sum("l", "q") + dims.SumSquared("s")

	if h.Rows() != cdim {
		err = errors.New(fmt.Sprintf("'h' must be float matrix of size (%d,1)", cdim))
		return
	}

	// Check b and set defaults if it is nil
	if b == nil {
		b = matrix.FloatZeros(0, 1)
	}
	if b.Cols() != 1 {
		estr := fmt.Sprintf("'b' must be a matrix with 1 column")
		err = errors.New(estr)
		return
	}

	mc := matrixVar{c}
	mb := matrixVar{b}
	var mG MatrixVarG
	var mA MatrixVarA

	if G == nil {
		mG = &matrixVarG{matrix.FloatZeros(0, c.Rows()), dims}
	} else {
		mG = &matrixIfG{G}
	}
	if A == nil {
		mA = &matrixVarA{matrix.FloatZeros(0, c.Rows())}
	} else {
		mA = &matrixIfA{A}
	}

	return cpl_problem(F, &mc, mG, h, mA, &mb, dims, kktsolver, solopts, x0, mnl)
}

func cpl_problem(F ConvexProg, c MatrixVariable, G MatrixVarG, h *matrix.FloatMatrix, A MatrixVarA,
	b MatrixVariable, dims *sets.DimensionSet, kktsolver KKTCpSolver,
	solopts *SolverOptions, x0 *matrix.FloatMatrix, mnl int) (sol *Solution, err error) {

	err = nil

	F_e := &convexVarProg{F}
	mx0 := &matrixVar{x0.Copy()}

	kktsolver_u := func(W *sets.FloatMatrixSet, x MatrixVariable, z *matrix.FloatMatrix) (KKTFuncVar, error) {
		g, err := kktsolver(W, x.Matrix(), z)
		solver := func(x, y MatrixVariable, z *matrix.FloatMatrix) error {
			return g(x.Matrix(), y.Matrix(), z)
		}
		return solver, err
	}
	return cpl_solver(F_e, c, G, h, A, b, dims, kktsolver_u, solopts, mx0, mnl)
}

// Internal CPL solver for CP and CLP problems. Everything is wrapped to proper interfaces
func cpl_solver(F ConvexVarProg, c MatrixVariable, G MatrixVarG, h *matrix.FloatMatrix,
	A MatrixVarA, b MatrixVariable, dims *sets.DimensionSet, kktsolver KKTCpSolverVar,
	solopts *SolverOptions, x0 MatrixVariable, mnl int) (sol *Solution, err error) {

	const (
		STEP              = 0.99
		BETA              = 0.5
		ALPHA             = 0.01
		EXPON             = 3
		MAX_RELAXED_ITERS = 8
	)

	var refinement int

	sol = &Solution{Unknown,
		nil,
		0.0, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 0}

	feasTolerance := FEASTOL
	absTolerance := ABSTOL
	relTolerance := RELTOL
	maxIter := MAXITERS
	if solopts.FeasTol > 0.0 {
		feasTolerance = solopts.FeasTol
	}
	if solopts.AbsTol > 0.0 {
		absTolerance = solopts.AbsTol
	}
	if solopts.RelTol > 0.0 {
		relTolerance = solopts.RelTol
	}
	if solopts.Refinement > 0 {
		refinement = solopts.Refinement
	} else {
		refinement = 1
	}
	if solopts.MaxIter > 0 {
		maxIter = solopts.MaxIter
	}

	if x0 == nil {
		mnl, x0, err = F.F0()
		if err != nil {
			return
		}
	}

	if c == nil {
		err = errors.New("Must define objective.")
		return
	}

	if h == nil {
		h = matrix.FloatZeros(0, 1)
	}
	if dims == nil {
		err = errors.New("Problem dimensions not defined.")
		return
	}
	if err = checkConeLpDimensions(dims); err != nil {
		return
	}

	cdim := dims.Sum("l", "q") + dims.SumSquared("s")
	cdim_diag := dims.Sum("l", "q", "s")

	if h.Rows() != cdim {
		err = errors.New(fmt.Sprintf("'h' must be float matrix of size (%d,1)", cdim))
		return
	}

	if G == nil {
		err = errors.New("'G' must be non-nil MatrixG interface.")
		return
	}
	fG := func(x, y MatrixVariable, alpha, beta float64, trans la.Option) error {
		return G.Gf(x, y, alpha, beta, trans)
	}

	// Check A and set defaults if it is nil
	if A == nil {
		err = errors.New("'A' must be non-nil MatrixA interface.")
		return
	}
	fA := func(x, y MatrixVariable, alpha, beta float64, trans la.Option) error {
		return A.Af(x, y, alpha, beta, trans)
	}

	if b == nil {
		err = errors.New("'b' must be non-nil MatrixVariable interface.")
		return
	}

	if kktsolver == nil {
		err = errors.New("nil kktsolver not allowed.")
		return
	}

	x := x0.Copy()
	y := b.Copy()
	y.Scal(0.0)
	z := matrix.FloatZeros(mnl+cdim, 1)
	s := matrix.FloatZeros(mnl+cdim, 1)
	ind := mnl + dims.At("l")[0]
	z.SetIndexes(1.0, matrix.MakeIndexSet(0, ind, 1)...)
	s.SetIndexes(1.0, matrix.MakeIndexSet(0, ind, 1)...)
	for _, m := range dims.At("q") {
		z.SetIndexes(1.0, ind)
		s.SetIndexes(1.0, ind)
		ind += m
	}
	for _, m := range dims.At("s") {
		iset := matrix.MakeIndexSet(ind, ind+m*m, m+1)
		z.SetIndexes(1.0, iset...)
		s.SetIndexes(1.0, iset...)
		ind += m * m
	}

	rx := x0.Copy()
	ry := b.Copy()
	dx := x.Copy()
	dy := y.Copy()
	rznl := matrix.FloatZeros(mnl, 1)
	rzl := matrix.FloatZeros(cdim, 1)
	dz := matrix.FloatZeros(mnl+cdim, 1)
	ds := matrix.FloatZeros(mnl+cdim, 1)
	lmbda := matrix.FloatZeros(mnl+cdim_diag, 1)
	lmbdasq := matrix.FloatZeros(mnl+cdim_diag, 1)
	sigs := matrix.FloatZeros(dims.Sum("s"), 1)
	sigz := matrix.FloatZeros(dims.Sum("s"), 1)

	dz2 := matrix.FloatZeros(mnl+cdim, 1)
	ds2 := matrix.FloatZeros(mnl+cdim, 1)

	newx := x.Copy()
	newy := y.Copy()
	newrx := x0.Copy()

	newz := matrix.FloatZeros(mnl+cdim, 1)
	news := matrix.FloatZeros(mnl+cdim, 1)
	newrznl := matrix.FloatZeros(mnl, 1)

	rx0 := rx.Copy()
	ry0 := ry.Copy()
	rznl0 := matrix.FloatZeros(mnl, 1)
	rzl0 := matrix.FloatZeros(cdim, 1)

	x0, dx0 := x.Copy(), dx.Copy()
	y0, dy0 := y.Copy(), dy.Copy()

	z0 := matrix.FloatZeros(mnl+cdim, 1)
	dz0 := matrix.FloatZeros(mnl+cdim, 1)
	dz20 := matrix.FloatZeros(mnl+cdim, 1)

	s0 := matrix.FloatZeros(mnl+cdim, 1)
	ds0 := matrix.FloatZeros(mnl+cdim, 1)
	ds20 := matrix.FloatZeros(mnl+cdim, 1)

	checkpnt.AddMatrixVar("z", z)
	checkpnt.AddMatrixVar("s", s)
	checkpnt.AddMatrixVar("dz", dz)
	checkpnt.AddMatrixVar("ds", ds)
	checkpnt.AddMatrixVar("rznl", rznl)
	checkpnt.AddMatrixVar("rzl", rzl)
	checkpnt.AddMatrixVar("lmbda", lmbda)
	checkpnt.AddMatrixVar("lmbdasq", lmbdasq)
	checkpnt.AddMatrixVar("z0", z0)
	checkpnt.AddMatrixVar("dz0", dz0)
	checkpnt.AddVerifiable("c", c)
	checkpnt.AddVerifiable("x", x)
	checkpnt.AddVerifiable("rx", rx)
	checkpnt.AddVerifiable("dx", dx)
	checkpnt.AddVerifiable("newrx", newrx)
	checkpnt.AddVerifiable("newx", newx)
	checkpnt.AddVerifiable("x0", x0)
	checkpnt.AddVerifiable("dx0", dx0)
	checkpnt.AddVerifiable("rx0", rx0)
	checkpnt.AddVerifiable("y", y)
	checkpnt.AddVerifiable("dy", dy)

	W0 := sets.NewFloatSet("d", "di", "dnl", "dnli", "v", "r", "rti", "beta")
	W0.Set("dnl", matrix.FloatZeros(mnl, 1))
	W0.Set("dnli", matrix.FloatZeros(mnl, 1))
	W0.Set("d", matrix.FloatZeros(dims.At("l")[0], 1))
	W0.Set("di", matrix.FloatZeros(dims.At("l")[0], 1))
	W0.Set("beta", matrix.FloatZeros(len(dims.At("q")), 1))
	for _, n := range dims.At("q") {
		W0.Append("v", matrix.FloatZeros(n, 1))
	}
	for _, n := range dims.At("s") {
		W0.Append("r", matrix.FloatZeros(n, n))
		W0.Append("rti", matrix.FloatZeros(n, n))
	}
	lmbda0 := matrix.FloatZeros(mnl+dims.Sum("l", "q", "s"), 1)
	lmbdasq0 := matrix.FloatZeros(mnl+dims.Sum("l", "q", "s"), 1)

	var f MatrixVariable = nil
	var Df MatrixVarDf = nil
	var H MatrixVarH = nil

	var ws3, wz3, wz2l, wz2nl *matrix.FloatMatrix
	var ws, wz, wz2, ws2 *matrix.FloatMatrix
	var wx, wx2, wy, wy2 MatrixVariable
	var gap, gap0, theta1, theta2, theta3, ts, tz, phi, phi0, mu, sigma, eta float64
	var resx, resy, reszl, resznl, pcost, dcost, dres, pres, relgap float64
	var resx0, resznl0, dres0, pres0 float64
	var dsdz, dsdz0, step, step0, dphi, dphi0, sigma0, eta0 float64
	var newresx, newresznl, newgap, newphi float64
	var W *sets.FloatMatrixSet
	var f3 KKTFuncVar

	checkpnt.AddFloatVar("gap", &gap)
	checkpnt.AddFloatVar("pcost", &pcost)
	checkpnt.AddFloatVar("dcost", &dcost)
	checkpnt.AddFloatVar("pres", &pres)
	checkpnt.AddFloatVar("dres", &dres)
	checkpnt.AddFloatVar("relgap", &relgap)
	checkpnt.AddFloatVar("step", &step)
	checkpnt.AddFloatVar("dsdz", &dsdz)
	checkpnt.AddFloatVar("resx", &resx)
	checkpnt.AddFloatVar("resy", &resy)
	checkpnt.AddFloatVar("reszl", &reszl)
	checkpnt.AddFloatVar("resznl", &resznl)

	// Declare fDf and fH here, they bind to Df and H as they are already declared.
	// ??really??

	var fDf func(u, v MatrixVariable, alpha, beta float64, trans la.Option) error = nil
	var fH func(u, v MatrixVariable, alpha, beta float64) error = nil

	relaxed_iters := 0
	for iters := 0; iters <= maxIter+1; iters++ {
		checkpnt.MajorNext()
		checkpnt.Check("loopstart", 10)

		checkpnt.MinorPush(10)
		if refinement != 0 || solopts.Debug {
			f, Df, H, err = F.F2(x, matrix.FloatVector(z.FloatArray()[:mnl]))
			fDf = func(u, v MatrixVariable, alpha, beta float64, trans la.Option) error {
				return Df.Df(u, v, alpha, beta, trans)
			}
			fH = func(u, v MatrixVariable, alpha, beta float64) error {
				return H.Hf(u, v, alpha, beta)
			}
		} else {
			f, Df, err = F.F1(x)
			fDf = func(u, v MatrixVariable, alpha, beta float64, trans la.Option) error {
				return Df.Df(u, v, alpha, beta, trans)
			}
		}
		checkpnt.MinorPop()

		gap = sdot(s, z, dims, mnl)

		// these are helpers, copies of parts of z,s
		z_mnl := matrix.FloatVector(z.FloatArray()[:mnl])
		z_mnl2 := matrix.FloatVector(z.FloatArray()[mnl:])
		s_mnl := matrix.FloatVector(s.FloatArray()[:mnl])
		s_mnl2 := matrix.FloatVector(s.FloatArray()[mnl:])

		// rx = c + A'*y + Df'*z[:mnl] + G'*z[mnl:]
		// -- y, rx MatrixArg
		mCopy(c, rx)
		fA(y, rx, 1.0, 1.0, la.OptTrans)
		fDf(&matrixVar{z_mnl}, rx, 1.0, 1.0, la.OptTrans)
		fG(&matrixVar{z_mnl2}, rx, 1.0, 1.0, la.OptTrans)
		resx = math.Sqrt(rx.Dot(rx))

		// rznl = s[:mnl] + f
		blas.Copy(s_mnl, rznl)
		blas.AxpyFloat(f.Matrix(), rznl, 1.0)
		resznl = blas.Nrm2Float(rznl)

		// rzl = s[mnl:] + G*x - h
		blas.Copy(s_mnl2, rzl)
		blas.AxpyFloat(h, rzl, -1.0)
		fG(x, &matrixVar{rzl}, 1.0, 1.0, la.OptNoTrans)
		reszl = snrm2(rzl, dims, 0)

		// Statistics for stopping criteria
		// pcost = c'*x
		// dcost = c'*x + y'*(A*x-b) + znl'*f(x) + zl'*(G*x-h)
		//       = c'*x + y'*(A*x-b) + znl'*(f(x)+snl) + zl'*(G*x-h+sl)
		//         - z'*s
		//       = c'*x + y'*ry + znl'*rznl + zl'*rzl - gap
		//pcost = blas.DotFloat(c, x)
		pcost = c.Dot(x)
		dcost = pcost + blas.DotFloat(y.Matrix(), ry.Matrix()) + blas.DotFloat(z_mnl, rznl)
		dcost += sdot(z_mnl2, rzl, dims, 0) - gap

		if pcost < 0.0 {
			relgap = gap / -pcost
		} else if dcost > 0.0 {
			relgap = gap / dcost
		} else {
			relgap = math.NaN()
		}
		pres = math.Sqrt(resy*resy + resznl*resznl + reszl*reszl)
		dres = resx
		if iters == 0 {
			resx0 = math.Max(1.0, resx)
			resznl0 = math.Max(1.0, resznl)
			pres0 = math.Max(1.0, pres)
			dres0 = math.Max(1.0, dres)
			gap0 = gap
			theta1 = 1.0 / gap0
			theta2 = 1.0 / resx0
			theta3 = 1.0 / resznl0
		}
		phi = theta1*gap + theta2*resx + theta3*resznl
		pres = pres / pres0
		dres = dres / dres0

		if solopts.ShowProgress {
			if iters == 0 {
				// some headers
				fmt.Printf("% 10s% 12s% 10s% 8s% 7s\n",
					"pcost", "dcost", "gap", "pres", "dres")
			}
			fmt.Printf("%2d: % 8.4e % 8.4e % 4.0e% 7.0e% 7.0e\n",
				iters, pcost, dcost, gap, pres, dres)
		}

		checkpnt.Check("checkgap", 50)
		// Stopping criteria
		if (pres <= feasTolerance && dres <= feasTolerance &&
			(gap <= absTolerance || (!math.IsNaN(relgap) && relgap <= relTolerance))) ||
			iters == maxIter {

			if iters == maxIter {
				s := "Terminated (maximum number of iterations reached)"
				if solopts.ShowProgress {
					fmt.Printf(s + "\n")
				}
				err = errors.New(s)
				sol.Status = Unknown
			} else {
				err = nil
				sol.Status = Optimal
			}
			sol.Result = sets.NewFloatSet("x", "y", "znl", "zl", "snl", "sl")
			sol.Result.Set("x", x.Matrix())
			sol.Result.Set("y", y.Matrix())
			sol.Result.Set("znl", matrix.FloatVector(z.FloatArray()[:mnl]))
			sol.Result.Set("zl", matrix.FloatVector(z.FloatArray()[mnl:]))
			sol.Result.Set("sl", matrix.FloatVector(s.FloatArray()[mnl:]))
			sol.Result.Set("snl", matrix.FloatVector(s.FloatArray()[:mnl]))
			sol.Gap = gap
			sol.RelativeGap = relgap
			sol.PrimalObjective = pcost
			sol.DualObjective = dcost
			sol.PrimalInfeasibility = pres
			sol.DualInfeasibility = dres
			sol.PrimalSlack = -ts
			sol.DualSlack = -tz
			return
		}

		// Compute initial scaling W:
		//
		//     W * z = W^{-T} * s = lambda.
		//
		// lmbdasq = lambda o lambda
		if iters == 0 {
			W, _ = computeScaling(s, z, lmbda, dims, mnl)
			checkpnt.AddScaleVar(W)
		}
		ssqr(lmbdasq, lmbda, dims, mnl)
		checkpnt.Check("lmbdasq", 90)

		// f3(x, y, z) solves
		//
		//     [ H   A'  GG'*W^{-1} ] [ ux ]   [ bx ]
		//     [ A   0   0          ] [ uy ] = [ by ].
		//     [ GG  0  -W'         ] [ uz ]   [ bz ]
		//
		// On entry, x, y, z contain bx, by, bz.
		// On exit, they contain ux, uy, uz.
		checkpnt.MinorPush(95)
		f3, err = kktsolver(W, x, z_mnl)
		checkpnt.MinorPop()
		checkpnt.Check("f3", 100)
		if err != nil {
			// ?? z_mnl is really copy of z[:mnl] ... should we copy here back to z??
			singular_kkt_matrix := false
			if iters == 0 {
				err = errors.New("Rank(A) < p or Rank([H(x); A; Df(x); G] < n")
				return
			} else if relaxed_iters > 0 && relaxed_iters < MAX_RELAXED_ITERS {
				// The arithmetic error may be caused by a relaxed line
				// search in the previous iteration.  Therefore we restore
				// the last saved state and require a standard line search.
				phi, gap = phi0, gap0
				mu = gap / float64(mnl+dims.Sum("l", "s")+len(dims.At("q")))
				blas.Copy(W0.At("dnl")[0], W.At("dnl")[0])
				blas.Copy(W0.At("dnli")[0], W.At("dnli")[0])
				blas.Copy(W0.At("d")[0], W.At("d")[0])
				blas.Copy(W0.At("di")[0], W.At("di")[0])
				blas.Copy(W0.At("beta")[0], W.At("beta")[0])
				for k, _ := range dims.At("q") {
					blas.Copy(W0.At("v")[k], W.At("v")[k])
				}
				for k, _ := range dims.At("s") {
					blas.Copy(W0.At("r")[k], W.At("r")[k])
					blas.Copy(W0.At("rti")[k], W.At("rti")[k])
				}
				//blas.Copy(x0, x)
				//x0.CopyTo(x)
				mCopy(x0, x)
				//blas.Copy(y0, y)
				mCopy(y0, y)
				blas.Copy(s0, s)
				blas.Copy(z0, z)
				blas.Copy(lmbda0, lmbda)
				blas.Copy(lmbdasq0, lmbdasq) // ???
				//blas.Copy(rx0, rx)
				//rx0.CopyTo(rx)
				mCopy(rx0, rx)
				//blas.Copy(ry0, ry)
				mCopy(ry0, ry)
				//resx = math.Sqrt(blas.DotFloat(rx, rx))
				resx = math.Sqrt(rx.Dot(rx))
				blas.Copy(rznl0, rznl)
				blas.Copy(rzl0, rzl)
				resznl = blas.Nrm2Float(rznl)

				relaxed_iters = -1

				// How about z_mnl here???
				checkpnt.MinorPush(120)
				f3, err = kktsolver(W, x, z_mnl)
				checkpnt.MinorPop()
				if err != nil {
					singular_kkt_matrix = true
				}
			} else {
				singular_kkt_matrix = true
			}

			if singular_kkt_matrix {
				msg := "Terminated (singular KKT matrix)."
				if solopts.ShowProgress {
					fmt.Printf(msg + "\n")
				}
				zl := matrix.FloatVector(z.FloatArray()[mnl:])
				sl := matrix.FloatVector(s.FloatArray()[mnl:])
				ind := dims.Sum("l", "q")
				for _, m := range dims.At("s") {
					symm(sl, m, ind)
					symm(zl, m, ind)
					ind += m * m
				}
				ts, _ = maxStep(s, dims, mnl, nil)
				tz, _ = maxStep(z, dims, mnl, nil)

				err = errors.New(msg)
				sol.Status = Unknown
				sol.Result = sets.NewFloatSet("x", "y", "znl", "zl", "snl", "sl")
				sol.Result.Set("x", x.Matrix())
				sol.Result.Set("y", y.Matrix())
				sol.Result.Set("znl", matrix.FloatVector(z.FloatArray()[:mnl]))
				sol.Result.Set("zl", zl)
				sol.Result.Set("sl", sl)
				sol.Result.Set("snl", matrix.FloatVector(s.FloatArray()[:mnl]))
				sol.Gap = gap
				sol.RelativeGap = relgap
				sol.PrimalObjective = pcost
				sol.DualObjective = dcost
				sol.PrimalInfeasibility = pres
				sol.DualInfeasibility = dres
				sol.PrimalSlack = -ts
				sol.DualSlack = -tz
				return
			}
		}

		// f4_no_ir(x, y, z, s) solves
		//
		//     [ 0     ]   [ H   A'  GG' ] [ ux        ]   [ bx ]
		//     [ 0     ] + [ A   0   0   ] [ uy        ] = [ by ]
		//     [ W'*us ]   [ GG  0   0   ] [ W^{-1}*uz ]   [ bz ]
		//
		//     lmbda o (uz + us) = bs.
		//
		// On entry, x, y, z, x, contain bx, by, bz, bs.
		// On exit, they contain ux, uy, uz, us.

		if iters == 0 {
			ws3 = matrix.FloatZeros(mnl+cdim, 1)
			wz3 = matrix.FloatZeros(mnl+cdim, 1)
			checkpnt.AddMatrixVar("ws3", ws3)
			checkpnt.AddMatrixVar("wz3", wz3)
		}

		f4_no_ir := func(x, y MatrixVariable, z, s *matrix.FloatMatrix) (err error) {
			// Solve
			//
			//     [ H  A'  GG'  ] [ ux        ]   [ bx                    ]
			//     [ A  0   0    ] [ uy        ] = [ by                    ]
			//     [ GG 0  -W'*W ] [ W^{-1}*uz ]   [ bz - W'*(lmbda o\ bs) ]
			//
			//     us = lmbda o\ bs - uz.

			err = nil
			// s := lmbda o\ s
			//    = lmbda o\ bs
			sinv(s, lmbda, dims, mnl)

			// z := z - W'*s
			//    = bz - W' * (lambda o\ bs)
			blas.Copy(s, ws3)

			scale(ws3, W, true, false)
			blas.AxpyFloat(ws3, z, -1.0)

			// Solve for ux, uy, uz
			err = f3(x, y, z)

			// s := s - z
			//    = lambda o\ bs - z.
			blas.AxpyFloat(z, s, -1.0)
			return
		}

		if iters == 0 {
			wz2nl = matrix.FloatZeros(mnl, 1)
			wz2l = matrix.FloatZeros(cdim, 1)
			checkpnt.AddMatrixVar("wz2nl", wz2nl)
			checkpnt.AddMatrixVar("wz2l", wz2l)
		}

		res := func(ux, uy MatrixVariable, uz, us *matrix.FloatMatrix, vx, vy MatrixVariable, vz, vs *matrix.FloatMatrix) (err error) {

			// Evaluates residuals in Newton equations:
			//
			//     [ vx ]     [ 0     ]   [ H  A' GG' ] [ ux        ]
			//     [ vy ] -=  [ 0     ] + [ A  0  0   ] [ uy        ]
			//     [ vz ]     [ W'*us ]   [ GG 0  0   ] [ W^{-1}*uz ]
			//
			//     vs -= lmbda o (uz + us).
			err = nil
			minor := checkpnt.MinorTop()
			// vx := vx - H*ux - A'*uy - GG'*W^{-1}*uz
			fH(ux, vx, -1.0, 1.0)
			fA(uy, vx, -1.0, 1.0, la.OptTrans)
			blas.Copy(uz, wz3)
			scale(wz3, W, false, true)
			wz3_nl := matrix.FloatVector(wz3.FloatArray()[:mnl])
			wz3_l := matrix.FloatVector(wz3.FloatArray()[mnl:])
			fDf(&matrixVar{wz3_nl}, vx, -1.0, 1.0, la.OptTrans)
			fG(&matrixVar{wz3_l}, vx, -1.0, 1.0, la.OptTrans)

			checkpnt.Check("10res", minor+10)

			// vy := vy - A*ux
			fA(ux, vy, -1.0, 1.0, la.OptNoTrans)

			// vz := vz - W'*us - GG*ux
			err = fDf(ux, &matrixVar{wz2nl}, 1.0, 0.0, la.OptNoTrans)
			checkpnt.Check("15res", minor+10)
			blas.AxpyFloat(wz2nl, vz, -1.0)
			fG(ux, &matrixVar{wz2l}, 1.0, 0.0, la.OptNoTrans)
			checkpnt.Check("20res", minor+10)
			blas.AxpyFloat(wz2l, vz, -1.0, &la.IOpt{"offsety", mnl})
			blas.Copy(us, ws3)
			scale(ws3, W, true, false)
			blas.AxpyFloat(ws3, vz, -1.0)

			checkpnt.Check("30res", minor+10)

			// vs -= lmbda o (uz + us)
			blas.Copy(us, ws3)
			blas.AxpyFloat(uz, ws3, 1.0)
			sprod(ws3, lmbda, dims, mnl, &la.SOpt{"diag", "D"})
			blas.AxpyFloat(ws3, vs, -1.0)

			checkpnt.Check("90res", minor+10)
			return
		}

		// f4(x, y, z, s) solves the same system as f4_no_ir, but applies
		// iterative refinement.

		if iters == 0 {
			if refinement > 0 || solopts.Debug {
				wx = c.Copy()
				wy = b.Copy()
				wz = z.Copy()
				ws = s.Copy()
				checkpnt.AddVerifiable("wx", wx)
				checkpnt.AddMatrixVar("ws", ws)
				checkpnt.AddMatrixVar("wz", wz)
			}
			if refinement > 0 {
				wx2 = c.Copy()
				wy2 = b.Copy()
				wz2 = matrix.FloatZeros(mnl+cdim, 1)
				ws2 = matrix.FloatZeros(mnl+cdim, 1)
				checkpnt.AddVerifiable("wx2", wx2)
				checkpnt.AddMatrixVar("ws2", ws2)
				checkpnt.AddMatrixVar("wz2", wz2)
			}
		}

		f4 := func(x, y MatrixVariable, z, s *matrix.FloatMatrix) (err error) {
			if refinement > 0 || solopts.Debug {
				mCopy(x, wx)
				mCopy(y, wy)
				blas.Copy(z, wz)
				blas.Copy(s, ws)
			}
			minor := checkpnt.MinorTop()
			checkpnt.Check("0_f4", minor+100)
			checkpnt.MinorPush(minor + 100)

			err = f4_no_ir(x, y, z, s)

			checkpnt.MinorPop()
			checkpnt.Check("1_f4", minor+200)
			for i := 0; i < refinement; i++ {
				mCopy(wx, wx2)
				mCopy(wy, wy2)
				blas.Copy(wz, wz2)
				blas.Copy(ws, ws2)

				checkpnt.Check("2_f4", minor+(1+i)*200)
				checkpnt.MinorPush(minor + (1+i)*200)

				res(x, y, z, s, wx2, wy2, wz2, ws2)
				checkpnt.MinorPop()
				checkpnt.Check("3_f4", minor+(1+i)*200+100)

				err = f4_no_ir(wx2, wy2, wz2, ws2)
				checkpnt.MinorPop()
				checkpnt.Check("4_f4", minor+(1+i)*200+199)
				wx2.Axpy(x, 1.0)
				wy2.Axpy(y, 1.0)
				blas.AxpyFloat(wz2, z, 1.0)
				blas.AxpyFloat(ws2, s, 1.0)
			}
			if solopts.Debug {
				res(x, y, z, s, wx, wy, wz, ws)
				fmt.Printf("KKT residuals:\n")
			}
			return
		}

		sigma, eta = 0.0, 0.0

		for i := 0; i < 2; i++ {
			minor := (i + 2) * 1000
			checkpnt.MinorPush(minor)
			checkpnt.Check("loop01", minor)

			// Solve
			//
			//     [ 0     ]   [ H  A' GG' ] [ dx        ]
			//     [ 0     ] + [ A  0  0   ] [ dy        ] = -(1 - eta)*r
			//     [ W'*ds ]   [ GG 0  0   ] [ W^{-1}*dz ]
			//
			//     lmbda o (dz + ds) = -lmbda o lmbda + sigma*mu*e.
			//

			mu = gap / float64(mnl+dims.Sum("l", "s")+len(dims.At("q")))
			blas.ScalFloat(ds, 0.0)
			blas.AxpyFloat(lmbdasq, ds, -1.0, &la.IOpt{"n", mnl + dims.Sum("l", "q")})

			ind = mnl + dims.At("l")[0]
			iset := matrix.MakeIndexSet(0, ind, 1)
			ds.Add(sigma*mu, iset...)
			for _, m := range dims.At("q") {
				ds.Add(sigma*mu, ind)
				ind += m
			}
			ind2 := ind
			for _, m := range dims.At("s") {
				blas.AxpyFloat(lmbdasq, ds, -1.0, &la.IOpt{"n", m}, &la.IOpt{"offsetx", ind2},
					&la.IOpt{"offsety", ind}, &la.IOpt{"incy", m + 1})
				ds.Add(sigma*mu, matrix.MakeIndexSet(ind, ind+m*m, m+1)...)
				ind += m * m
				ind2 += m
			}

			dx.Scal(0.0)
			rx.Axpy(dx, -1.0+eta)
			dy.Scal(0.0)
			ry.Axpy(dy, -1.0+eta)
			dz.Scale(0.0)
			blas.AxpyFloat(rznl, dz, -1.0+eta)
			blas.AxpyFloat(rzl, dz, -1.0+eta, &la.IOpt{"offsety", mnl})
			//fmt.Printf("dx=\n%v\n", dx)
			//fmt.Printf("dz=\n%v\n", dz.ToString("%.7f"))
			//fmt.Printf("ds=\n%v\n", ds.ToString("%.7f"))

			checkpnt.Check("pref4", minor)
			checkpnt.MinorPush(minor)
			err = f4(dx, dy, dz, ds)
			if err != nil {
				if iters == 0 {
					s := fmt.Sprintf("Rank(A) < p or Rank([H(x); A; Df(x); G] < n (%s)", err)
					err = errors.New(s)
					return
				}
				msg := "Terminated (singular KKT matrix)."
				if solopts.ShowProgress {
					fmt.Printf(msg + "\n")
				}
				zl := matrix.FloatVector(z.FloatArray()[mnl:])
				sl := matrix.FloatVector(s.FloatArray()[mnl:])
				ind := dims.Sum("l", "q")
				for _, m := range dims.At("s") {
					symm(sl, m, ind)
					symm(zl, m, ind)
					ind += m * m
				}
				ts, _ = maxStep(s, dims, mnl, nil)
				tz, _ = maxStep(z, dims, mnl, nil)

				err = errors.New(msg)
				sol.Status = Unknown
				sol.Result = sets.NewFloatSet("x", "y", "znl", "zl", "snl", "sl")
				sol.Result.Set("x", x.Matrix())
				sol.Result.Set("y", y.Matrix())
				sol.Result.Set("znl", matrix.FloatVector(z.FloatArray()[:mnl]))
				sol.Result.Set("zl", zl)
				sol.Result.Set("sl", sl)
				sol.Result.Set("snl", matrix.FloatVector(s.FloatArray()[:mnl]))
				sol.Gap = gap
				sol.RelativeGap = relgap
				sol.PrimalObjective = pcost
				sol.DualObjective = dcost
				sol.PrimalInfeasibility = pres
				sol.DualInfeasibility = dres
				sol.PrimalSlack = -ts
				sol.DualSlack = -tz
				return
			}

			checkpnt.MinorPop()
			checkpnt.Check("postf4", minor+400)

			// Inner product ds'*dz and unscaled steps are needed in the
			// line search.
			dsdz = sdot(ds, dz, dims, mnl)
			blas.Copy(dz, dz2)
			scale(dz2, W, false, true)
			blas.Copy(ds, ds2)
			scale(ds2, W, true, false)

			checkpnt.Check("dsdz", minor+400)

			// Maximum steps to boundary.
			//
			// Also compute the eigenvalue decomposition of 's' blocks in
			// ds, dz.  The eigenvectors Qs, Qz are stored in ds, dz.
			// The eigenvalues are stored in sigs, sigz.

			scale2(lmbda, ds, dims, mnl, false)
			ts, _ = maxStep(ds, dims, mnl, sigs)
			scale2(lmbda, dz, dims, mnl, false)
			tz, _ = maxStep(dz, dims, mnl, sigz)
			t := maxvec([]float64{0.0, ts, tz})
			if t == 0 {
				step = 1.0
			} else {
				step = math.Min(1.0, STEP/t)
			}

			checkpnt.Check("maxstep", minor+400)

			var newDf MatrixVarDf = nil
			var newf MatrixVariable = nil

			// Backtrack until newx is in domain of f.
			backtrack := true
			for backtrack {
				mCopy(x, newx)
				dx.Axpy(newx, step)
				newf, newDf, err = F.F1(newx)
				if newf != nil {
					backtrack = false
				} else {
					step *= BETA
				}
			}

			// Merit function
			//
			//     phi = theta1 * gap + theta2 * norm(rx) +
			//         theta3 * norm(rznl)
			//
			// and its directional derivative dphi.

			phi = theta1*gap + theta2*resx + theta3*resznl
			if i == 0 {
				dphi = -phi
			} else {
				dphi = -theta1*(1-sigma)*gap - theta2*(1-eta)*resx - theta3*(1-eta)*resznl
			}

			var newfDf func(x, y MatrixVariable, a, b float64, trans la.Option) error

			// Line search
			backtrack = true
			for backtrack {
				mCopy(x, newx)
				dx.Axpy(newx, step)
				mCopy(y, newy)
				dy.Axpy(newy, step)
				blas.Copy(z, newz)
				blas.AxpyFloat(dz2, newz, step)
				blas.Copy(s, news)
				blas.AxpyFloat(ds2, news, step)

				newf, newDf, err = F.F1(newx)
				newfDf = func(u, v MatrixVariable, a, b float64, trans la.Option) error {
					return newDf.Df(u, v, a, b, trans)
				}

				// newrx = c + A'*newy + newDf'*newz[:mnl] + G'*newz[mnl:]
				newz_mnl := matrix.FloatVector(newz.FloatArray()[:mnl])
				newz_ml := matrix.FloatVector(newz.FloatArray()[mnl:])
				//blas.Copy(c, newrx)
				//c.CopyTo(newrx)
				mCopy(c, newrx)
				fA(newy, newrx, 1.0, 1.0, la.OptTrans)
				newfDf(&matrixVar{newz_mnl}, newrx, 1.0, 1.0, la.OptTrans)
				fG(&matrixVar{newz_ml}, newrx, 1.0, 1.0, la.OptTrans)
				newresx = math.Sqrt(newrx.Dot(newrx))

				// newrznl = news[:mnl] + newf
				news_mnl := matrix.FloatVector(news.FloatArray()[:mnl])
				//news_ml := matrix.FloatVector(news.FloatArray()[mnl:])
				blas.Copy(news_mnl, newrznl)
				blas.AxpyFloat(newf.Matrix(), newrznl, 1.0)
				newresznl = blas.Nrm2Float(newrznl)

				newgap = (1.0-(1.0-sigma)*step)*gap + step*step*dsdz
				newphi = theta1*newgap + theta2*newresx + theta3*newresznl

				if i == 0 {
					if newgap <= (1.0-ALPHA*step)*gap &&
						(relaxed_iters > 0 && relaxed_iters < MAX_RELAXED_ITERS ||
							newphi <= phi+ALPHA*step*dphi) {
						backtrack = false
						sigma = math.Min(newgap/gap, math.Pow((newgap/gap), EXPON))
						//fmt.Printf("break 1: sigma=%.7f\n", sigma)
						eta = 0.0
					} else {
						step *= BETA
					}
				} else {
					if relaxed_iters == -1 || (relaxed_iters == 0 && MAX_RELAXED_ITERS == 0) {
						// Do a standard line search.
						if newphi <= phi+ALPHA*step*dphi {
							relaxed_iters = 0
							backtrack = false
							//fmt.Printf("break 2 : newphi=%.7f\n", newphi)
						} else {
							step *= BETA
						}
					} else if relaxed_iters == 0 && relaxed_iters < MAX_RELAXED_ITERS {
						if newphi <= phi+ALPHA*step*dphi {
							// Relaxed l.s. gives sufficient decrease.
							relaxed_iters = 0
						} else {
							// Save state.
							phi0, dphi0, gap0 = phi, dphi, gap
							step0 = step

							blas.Copy(W.At("dnl")[0], W0.At("dnl")[0])
							blas.Copy(W.At("dnli")[0], W0.At("dnli")[0])
							blas.Copy(W.At("d")[0], W0.At("d")[0])
							blas.Copy(W.At("di")[0], W0.At("di")[0])
							blas.Copy(W.At("beta")[0], W0.At("beta")[0])
							for k, _ := range dims.At("q") {
								blas.Copy(W.At("v")[k], W0.At("v")[k])
							}
							for k, _ := range dims.At("s") {
								blas.Copy(W.At("r")[k], W0.At("r")[k])
								blas.Copy(W.At("rti")[k], W0.At("rti")[k])
							}
							mCopy(x, x0)
							mCopy(y, y0)
							mCopy(dx, dx0)
							mCopy(dy, dy0)
							blas.Copy(s, s0)
							blas.Copy(z, z0)
							blas.Copy(ds, ds0)
							blas.Copy(dz, dz0)
							blas.Copy(ds2, ds20)
							blas.Copy(dz2, dz20)
							blas.Copy(lmbda, lmbda0)
							blas.Copy(lmbdasq, lmbdasq0) // ???
							mCopy(rx, rx0)
							mCopy(ry, ry0)
							blas.Copy(rznl, rznl0)
							blas.Copy(rzl, rzl0)
							dsdz0 = dsdz
							sigma0, eta0 = sigma, eta
							relaxed_iters = 1
						}
						backtrack = false
						//fmt.Printf("break 3 : newphi=%.7f\n", newphi)

					} else if relaxed_iters >= 0 && relaxed_iters < MAX_RELAXED_ITERS &&
						MAX_RELAXED_ITERS > 0 {
						if newphi <= phi0+ALPHA*step0*dphi0 {
							// Relaxed l.s. gives sufficient decrease.
							relaxed_iters = 0
						} else {
							// Relaxed line search
							relaxed_iters += 1
						}
						backtrack = false
						//fmt.Printf("break 4 : newphi=%.7f\n", newphi)

					} else if relaxed_iters == MAX_RELAXED_ITERS && MAX_RELAXED_ITERS > 0 {
						if newphi <= phi0+ALPHA*step0*dphi0 {
							// Series of relaxed line searches ends
							// with sufficient decrease w.r.t. phi0.
							backtrack = false
							relaxed_iters = 0
							//fmt.Printf("break 5 : newphi=%.7f\n", newphi)
						} else if newphi >= phi0 {
							// Resume last saved line search
							phi, dphi, gap = phi0, dphi0, gap0
							step = step0
							blas.Copy(W0.At("dnl")[0], W.At("dnl")[0])
							blas.Copy(W0.At("dnli")[0], W.At("dnli")[0])
							blas.Copy(W0.At("d")[0], W.At("d")[0])
							blas.Copy(W0.At("di")[0], W.At("di")[0])
							blas.Copy(W0.At("beta")[0], W.At("beta")[0])
							for k, _ := range dims.At("q") {
								blas.Copy(W0.At("v")[k], W.At("v")[k])
							}
							for k, _ := range dims.At("s") {
								blas.Copy(W0.At("r")[k], W.At("r")[k])
								blas.Copy(W0.At("rti")[k], W.At("rti")[k])
							}
							mCopy(x, x0)
							mCopy(y, y0)
							mCopy(dx, dx0)
							mCopy(dy, dy0)
							blas.Copy(s, s0)
							blas.Copy(z, z0)
							blas.Copy(ds2, ds20)
							blas.Copy(dz2, dz20)
							blas.Copy(lmbda, lmbda0)
							blas.Copy(lmbdasq, lmbdasq0) // ???
							mCopy(rx, rx0)
							mCopy(ry, ry0)
							blas.Copy(rznl, rznl0)
							blas.Copy(rzl, rzl0)
							dsdz = dsdz0
							sigma, eta = sigma0, eta0
							relaxed_iters = -1

						} else if newphi <= phi+ALPHA*step*dphi {
							// Series of relaxed line searches ends
							// with sufficient decrease w.r.t. phi0.
							backtrack = false
							relaxed_iters = -1
							//fmt.Printf("break 6 : newphi=%.7f\n", newphi)
						}
					}
				}
			} // end of line search

			checkpnt.Check("eol", minor+900)

		} // end for [0,1]

		// Update x, y
		dx.Axpy(x, step)
		dy.Axpy(y, step)
		checkpnt.Check("updatexy", 5000)

		// Replace nonlinear, 'l' and 'q' blocks of ds and dz with the
		// updated variables in the current scaling.
		// Replace 's' blocks of ds and dz with the factors Ls, Lz in a
		// factorization Ls*Ls', Lz*Lz' of the updated variables in the
		// current scaling.

		// ds := e + step*ds for nonlinear, 'l' and 'q' blocks.
		// dz := e + step*dz for nonlinear, 'l' and 'q' blocks.
		blas.ScalFloat(ds, step, &la.IOpt{"n", mnl + dims.Sum("l", "q")})
		blas.ScalFloat(dz, step, &la.IOpt{"n", mnl + dims.Sum("l", "q")})
		ind := mnl + dims.At("l")[0]
		is := matrix.MakeIndexSet(0, ind, 1)
		ds.Add(1.0, is...)
		dz.Add(1.0, is...)
		for _, m := range dims.At("q") {
			ds.SetIndex(ind, 1.0+ds.GetIndex(ind))
			dz.SetIndex(ind, 1.0+dz.GetIndex(ind))
			ind += m
		}
		checkpnt.Check("updatedsdz", 5100)

		// ds := H(lambda)^{-1/2} * ds and dz := H(lambda)^{-1/2} * dz.
		//
		// This replaces the 'l' and 'q' components of ds and dz with the
		// updated variables in the current scaling.
		// The 's' components of ds and dz are replaced with
		//
		// diag(lmbda_k)^{1/2} * Qs * diag(lmbda_k)^{1/2}
		// diag(lmbda_k)^{1/2} * Qz * diag(lmbda_k)^{1/2}
		scale2(lmbda, ds, dims, mnl, true)
		scale2(lmbda, dz, dims, mnl, true)

		checkpnt.Check("scale2", 5200)

		// sigs := ( e + step*sigs ) ./ lambda for 's' blocks.
		// sigz := ( e + step*sigz ) ./ lambda for 's' blocks.
		blas.ScalFloat(sigs, step)
		blas.ScalFloat(sigz, step)
		sigs.Add(1.0)
		sigz.Add(1.0)
		sdimsum := dims.Sum("s")
		qdimsum := dims.Sum("l", "q")
		blas.TbsvFloat(lmbda, sigs, &la.IOpt{"n", sdimsum}, &la.IOpt{"k", 0},
			&la.IOpt{"lda", 1}, &la.IOpt{"offseta", mnl + qdimsum})
		blas.TbsvFloat(lmbda, sigz, &la.IOpt{"n", sdimsum}, &la.IOpt{"k", 0},
			&la.IOpt{"lda", 1}, &la.IOpt{"offseta", mnl + qdimsum})

		checkpnt.Check("sigs", 5300)

		ind2 := mnl + qdimsum
		ind3 := 0
		sdims := dims.At("s")

		for k := 0; k < len(sdims); k++ {
			m := sdims[k]
			for i := 0; i < m; i++ {
				a := math.Sqrt(sigs.GetIndex(ind3 + i))
				blas.ScalFloat(ds, a, &la.IOpt{"offset", ind2 + m*i}, &la.IOpt{"n", m})
				a = math.Sqrt(sigz.GetIndex(ind3 + i))
				blas.ScalFloat(dz, a, &la.IOpt{"offset", ind2 + m*i}, &la.IOpt{"n", m})
			}
			ind2 += m * m
			ind3 += m
		}

		checkpnt.Check("scaling", 5400)
		err = updateScaling(W, lmbda, ds, dz)
		checkpnt.Check("postscaling", 5500)

		// Unscale s, z, tau, kappa (unscaled variables are used only to
		// compute feasibility residuals).
		ind = mnl + dims.Sum("l", "q")
		ind2 = ind
		blas.Copy(lmbda, s, &la.IOpt{"n", ind})
		for _, m := range dims.At("s") {
			blas.ScalFloat(s, 0.0, &la.IOpt{"offset", ind2})
			blas.Copy(lmbda, s, &la.IOpt{"offsetx", ind}, &la.IOpt{"offsety", ind2},
				&la.IOpt{"n", m}, &la.IOpt{"incy", m + 1})
			ind += m
			ind2 += m * m
		}
		scale(s, W, true, false)
		checkpnt.Check("unscale_s", 5600)

		ind = mnl + dims.Sum("l", "q")
		ind2 = ind
		blas.Copy(lmbda, z, &la.IOpt{"n", ind})
		for _, m := range dims.At("s") {
			blas.ScalFloat(z, 0.0, &la.IOpt{"offset", ind2})
			blas.Copy(lmbda, z, &la.IOpt{"offsetx", ind}, &la.IOpt{"offsety", ind2},
				&la.IOpt{"n", m}, &la.IOpt{"incy", m + 1})
			ind += m
			ind2 += m * m
		}
		scale(z, W, false, true)
		checkpnt.Check("unscale_z", 5700)

		gap = blas.DotFloat(lmbda, lmbda)

	}
	return
}

// Local Variables:
// tab-width: 4
// End:
