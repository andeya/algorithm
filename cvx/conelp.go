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

// Implements MatrixA interface for standard matrix valued A.
type matrixA struct {
	mA *matrix.FloatMatrix
}

func (a *matrixA) Af(x, y *matrix.FloatMatrix, alpha, beta float64, trans la.Option) error {
	return blas.GemvFloat(a.mA, x, y, alpha, beta, trans)
}

// Implements MatrixG interface for standard matrix valued G.
type matrixG struct {
	mG   *matrix.FloatMatrix
	dims *sets.DimensionSet
}

func (g *matrixG) Gf(x, y *matrix.FloatMatrix, alpha, beta float64, trans la.Option) error {
	return sgemv(g.mG, x, y, alpha, beta, g.dims, trans)
}

type fClosure struct {
	wx, wy, ws, wz     *matrix.FloatMatrix
	wx2, wy2, ws2, wz2 *matrix.FloatMatrix
	// these are singleton matrices
	wtau, wkappa, wtau2, wkappa2 *matrix.FloatMatrix
}

type fVarClosure struct {
	wx, wy   MatrixVariable
	ws, wz   *matrix.FloatMatrix
	wx2, wy2 MatrixVariable
	ws2, wz2 *matrix.FloatMatrix
	// these are singleton matrices
	wtau, wkappa, wtau2, wkappa2 *matrix.FloatMatrix
}

func checkConeLpDimensions(dims *sets.DimensionSet) error {
	if len(dims.At("l")) == 0 {
		dims.Set("l", []int{0})
	} else if dims.At("l")[0] < 0 {
		return errors.New("dimension 'l' must be nonnegative integer")
	}
	for _, m := range dims.At("q") {
		if m < 1 {
			return errors.New("dimension 'q' must be list of positive integers")
		}
	}
	for _, m := range dims.At("s") {
		if m < 1 {
			return errors.New("dimension 's' must be list of positive integers")
		}
	}
	return nil
}

// Solves a pair of primal and dual cone programs
//
//        minimize    c'*x
//        subject to  G*x + s = h
//                    A*x = b
//                    s >= 0
//
//        maximize    -h'*z - b'*y
//        subject to  G'*z + A'*y + c = 0
//                    z >= 0.
//
// The inequalities are with respect to a cone C defined as the Cartesian
// product of N + M + 1 cones:
//
//        C = C_0 x C_1 x .... x C_N x C_{N+1} x ... x C_{N+M}.
//
// The first cone C_0 is the nonnegative orthant of dimension ml.
// The next N cones are second order cones of dimension r[0], ..., r[N-1].
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
// The default value for dims is l: []int{G.Rows()}, q: []int{}, s: []int{}.
//
// Arguments primalstart, dualstart are optional starting points for primal and
// dual problems. If non-nil then primalstart is a FloatMatrixSet having two entries.
//
//  primalstart.At("x")[0]  starting point for x
//  primalstart.At("s")[0]  starting point for s
//  dualstart.At("y")[0]    starting point for y
//  dualstart.At("z")[0]    starting point for z
//
// On exit Solution contains the result and information about the accurancy of the
// solution. if SolutionStatus is Optimal then Solution.Result contains solutions
// for the problems.
//
//   Result.At("x")[0]  solution for x
//   Result.At("y")[0]  solution for y
//   Result.At("s")[0]  solution for s
//   Result.At("z")[0]  solution for z
//
func ConeLp(c, G, h, A, b *matrix.FloatMatrix, dims *sets.DimensionSet, solopts *SolverOptions,
	primalstart, dualstart *sets.FloatMatrixSet) (sol *Solution, err error) {

	if c == nil || c.Cols() > 1 {
		err = errors.New("'c' must be matrix with 1 column")
		return
	}
	if c.Rows() < 1 {
		err = errors.New("No variables, 'c' must have at least one row")
		return

	}
	if h == nil || h.Cols() > 1 {
		err = errors.New("'h' must be matrix with 1 column")
		return
	}

	if dims == nil {
		dims = sets.NewDimensionSet("l", "q", "s")
		dims.Set("l", []int{h.Rows()})
	}

	cdim := dims.Sum("l", "q") + dims.SumSquared("s")
	cdim_pckd := dims.Sum("l", "q") + dims.SumPacked("s")

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

	if b.Rows() > c.Rows() || b.Rows()+cdim_pckd < c.Rows() {
		err = errors.New("Rank(A) < p or Rank([G; A]) < n")
		return
	}

	solvername := solopts.KKTSolverName
	if len(solvername) == 0 {
		if len(dims.At("q")) > 0 || len(dims.At("s")) > 0 {
			solvername = "qr"
		} else {
			solvername = "chol2"
		}
	}

	var factor kktFactor
	var kktsolver KKTConeSolver = nil
	if kktfunc, ok := lpsolvers[solvername]; ok {
		// kkt function returns us problem spesific factor function.
		factor, err = kktfunc(G, dims, A, 0)
		if err != nil {
			return nil, err
		}
		kktsolver = func(W *sets.FloatMatrixSet) (KKTFunc, error) {
			return factor(W, nil, nil)
		}
	} else {
		err = errors.New(fmt.Sprintf("solver '%s' not known", solvername))
		return
	}
	//return ConeLpCustom(c, &mG, h, &mA, b, dims, kktsolver, solopts, primalstart, dualstart)
	c_e := &matrixVar{c}
	G_e := &matrixVarG{G, dims}
	A_e := &matrixVarA{A}
	b_e := &matrixVar{b}
	return conelp_problem(c_e, G_e, h, A_e, b_e, dims, kktsolver, solopts, primalstart, dualstart)
}

// Solves a pair of primal and dual cone programs  using custom KKT solver.
//
func ConeLpCustomKKT(c, G, h, A, b *matrix.FloatMatrix, dims *sets.DimensionSet,
	kktsolver KKTConeSolver, solopts *SolverOptions, primalstart,
	dualstart *sets.FloatMatrixSet) (sol *Solution, err error) {

	if c == nil || c.Cols() > 1 {
		err = errors.New("'c' must be matrix with 1 column")
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
	cdim_pckd := dims.Sum("l", "q") + dims.SumPacked("s")
	//cdim_diag := dims.Sum("l", "q", "s")

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

	if b.Rows() > c.Rows() || b.Rows()+cdim_pckd < c.Rows() {
		err = errors.New("Rank(A) < p or Rank([G; A]) < n")
		return
	}

	mA := &matrixVarA{A}
	mG := &matrixVarG{G, dims}
	mc := &matrixVar{c}
	mb := &matrixVar{b}

	return conelp_problem(mc, mG, h, mA, mb, dims, kktsolver, solopts, primalstart, dualstart)
}

// Solves a pair of primal and dual cone programs using custom KKT solver and constraint
// interfaces MatrixG and MatrixA
//
func ConeLpCustomMatrix(c *matrix.FloatMatrix, G MatrixG, h *matrix.FloatMatrix,
	A MatrixA, b *matrix.FloatMatrix, dims *sets.DimensionSet, kktsolver KKTConeSolver,
	solopts *SolverOptions, primalstart, dualstart *sets.FloatMatrixSet) (sol *Solution, err error) {

	err = nil

	if c == nil || c.Cols() > 1 {
		err = errors.New("'c' must be matrix with 1 column")
		return
	}
	if h == nil || h.Cols() > 1 {
		err = errors.New("'h' must be matrix with 1 column")
		return
	}

	if err = checkConeLpDimensions(dims); err != nil {
		return
	}

	cdim := dims.Sum("l", "q") + dims.SumSquared("s")
	cdim_pckd := dims.Sum("l", "q") + dims.SumPacked("s")
	//cdim_diag := dims.Sum("l", "q", "s")

	if h.Rows() != cdim {
		err = errors.New(fmt.Sprintf("'h' must be float matrix of size (%d,1)", cdim))
		return
	}

	// Data for kth 'q' constraint are found in rows indq[k]:indq[k+1] of G.
	indq := make([]int, 0)
	indq = append(indq, dims.At("l")[0])
	for _, k := range dims.At("q") {
		indq = append(indq, indq[len(indq)-1]+k)
	}

	// Data for kth 's' constraint are found in rows inds[k]:inds[k+1] of G.
	inds := make([]int, 0)
	inds = append(inds, indq[len(indq)-1])
	for _, k := range dims.At("s") {
		inds = append(inds, inds[len(inds)-1]+k*k)
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
	if b.Rows() > c.Rows() || b.Rows()+cdim_pckd < c.Rows() {
		err = errors.New("Rank(A) < p or Rank([G; A]) < n")
		return
	}

	if kktsolver == nil {
		err = errors.New("nil kktsolver not allowed.")
		return
	}

	var mA MatrixVarA
	var mG MatrixVarG
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
	var mc = &matrixVar{c}
	var mb = &matrixVar{b}

	return conelp_problem(mc, mG, h, mA, mb, dims, kktsolver, solopts, primalstart, dualstart)
}

func conelp_problem(c MatrixVariable, G MatrixVarG, h *matrix.FloatMatrix,
	A MatrixVarA, b MatrixVariable, dims *sets.DimensionSet, kktsolver KKTConeSolver,
	solopts *SolverOptions, primalstart, dualstart *sets.FloatMatrixSet) (sol *Solution, err error) {

	kktsolver_u := func(W *sets.FloatMatrixSet) (KKTFuncVar, error) {
		g, err := kktsolver(W)
		solver := func(x, y MatrixVariable, z *matrix.FloatMatrix) error {
			return g(x.Matrix(), y.Matrix(), z)
		}
		return solver, err
	}
	return conelp_solver(c, G, h, A, b, dims, kktsolver_u, solopts, primalstart, dualstart)
}

func conelp_solver(c MatrixVariable, G MatrixVarG, h *matrix.FloatMatrix,
	A MatrixVarA, b MatrixVariable, dims *sets.DimensionSet, kktsolver KKTConeSolverVar,
	solopts *SolverOptions, primalstart, dualstart *sets.FloatMatrixSet) (sol *Solution, err error) {

	err = nil
	const EXPON = 3
	const STEP = 0.99

	sol = &Solution{Unknown,
		nil,
		0.0, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 0}

	var refinement int

	if solopts.Refinement > 0 {
		refinement = solopts.Refinement
	} else {
		refinement = 0
		if len(dims.At("q")) > 0 || len(dims.At("s")) > 0 {
			refinement = 1
		}
	}
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
	if solopts.MaxIter > 0 {
		maxIter = solopts.MaxIter
	}
	if err = checkConeLpDimensions(dims); err != nil {
		return
	}

	cdim := dims.Sum("l", "q") + dims.SumSquared("s")
	//cdim_pckd := dims.Sum("l", "q") + dims.SumPacked("s")
	cdim_diag := dims.Sum("l", "q", "s")

	if h.Rows() != cdim {
		err = errors.New(fmt.Sprintf("'h' must be float matrix of size (%d,1)", cdim))
		return
	}

	// Data for kth 'q' constraint are found in rows indq[k]:indq[k+1] of G.
	indq := make([]int, 0)
	indq = append(indq, dims.At("l")[0])
	for _, k := range dims.At("q") {
		indq = append(indq, indq[len(indq)-1]+k)
	}

	// Data for kth 's' constraint are found in rows inds[k]:inds[k+1] of G.
	inds := make([]int, 0)
	inds = append(inds, indq[len(indq)-1])
	for _, k := range dims.At("s") {
		inds = append(inds, inds[len(inds)-1]+k*k)
	}

	Gf := func(x, y MatrixVariable, alpha, beta float64, trans la.Option) error {
		return G.Gf(x, y, alpha, beta, trans)
	}

	Af := func(x, y MatrixVariable, alpha, beta float64, trans la.Option) error {
		return A.Af(x, y, alpha, beta, trans)
	}

	// kktsolver(W) returns a routine for solving 3x3 block KKT system
	//
	//     [ 0   A'  G'*W^{-1} ] [ ux ]   [ bx ]
	//     [ A   0   0         ] [ uy ] = [ by ].
	//     [ G   0   -W'       ] [ uz ]   [ bz ]

	if kktsolver == nil {
		err = errors.New("nil kktsolver not allowed.")
		return
	}

	// res() evaluates residual in 5x5 block KKT system
	//
	//     [ vx   ]    [ 0         ]   [ 0   A'  G'  c ] [ ux        ]
	//     [ vy   ]    [ 0         ]   [-A   0   0   b ] [ uy        ]
	//     [ vz   ] += [ W'*us     ] - [-G   0   0   h ] [ W^{-1}*uz ]
	//     [ vtau ]    [ dg*ukappa ]   [-c' -b' -h'  0 ] [ utau/dg   ]
	//
	//           vs += lmbda o (dz + ds)
	//       vkappa += lmbdg * (dtau + dkappa).
	ws3 := matrix.FloatZeros(cdim, 1)
	wz3 := matrix.FloatZeros(cdim, 1)
	checkpnt.AddMatrixVar("ws3", ws3)
	checkpnt.AddMatrixVar("wz3", wz3)

	//
	res := func(ux, uy MatrixVariable, uz, utau, us, ukappa *matrix.FloatMatrix,
		vx, vy MatrixVariable, vz, vtau, vs, vkappa *matrix.FloatMatrix,
		W *sets.FloatMatrixSet, dg float64, lmbda *matrix.FloatMatrix) (err error) {

		err = nil
		// vx := vx - A'*uy - G'*W^{-1}*uz - c*utau/dg
		Af(uy, vx, -1.0, 1.0, la.OptTrans)
		//fmt.Printf("post-Af vx=\n%v\n", vx)
		blas.Copy(uz, wz3)
		scale(wz3, W, false, true)
		Gf(&matrixVar{wz3}, vx, -1.0, 1.0, la.OptTrans)
		//blas.AxpyFloat(c, vx, -utau.Float()/dg)
		c.Axpy(vx, -utau.Float()/dg)

		// vy := vy + A*ux - b*utau/dg
		Af(ux, vy, 1.0, 1.0, la.OptNoTrans)
		//blas.AxpyFloat(b, vy, -utau.Float()/dg)
		b.Axpy(vy, -utau.Float()/dg)

		// vz := vz + G*ux - h*utau/dg + W'*us
		Gf(ux, &matrixVar{vz}, 1.0, 1.0, la.OptNoTrans)
		blas.AxpyFloat(h, vz, -utau.Float()/dg)
		blas.Copy(us, ws3)
		scale(ws3, W, true, false)
		blas.AxpyFloat(ws3, vz, 1.0)

		// vtau := vtau + c'*ux + b'*uy + h'*W^{-1}*uz + dg*ukappa
		var vtauplus float64 = dg*ukappa.Float() + c.Dot(ux) +
			b.Dot(uy) + sdot(h, wz3, dims, 0)
		vtau.SetValue(vtau.Float() + vtauplus)

		// vs := vs + lmbda o (uz + us)
		blas.Copy(us, ws3)
		blas.AxpyFloat(uz, ws3, 1.0)
		sprod(ws3, lmbda, dims, 0, &la.SOpt{"diag", "D"})
		blas.AxpyFloat(ws3, vs, 1.0)

		// vkappa += vkappa + lmbdag * (utau + ukappa)
		lscale := lmbda.GetIndex(-1)
		var vkplus float64 = lscale * (utau.Float() + ukappa.Float())
		vkappa.SetValue(vkappa.Float() + vkplus)
		return
	}

	resx0 := math.Max(1.0, math.Sqrt(c.Dot(c)))
	resy0 := math.Max(1.0, math.Sqrt(b.Dot(b)))
	resz0 := math.Max(1.0, snrm2(h, dims, 0))

	// select initial points

	//fmt.Printf("** initial resx0=%.4f, resy0=%.4f, resz0=%.4f \n", resx0, resy0, resz0)

	x := c.Copy()
	//blas.ScalFloat(x, 0.0)
	x.Scal(0.0)
	y := b.Copy()
	//blas.ScalFloat(y, 0.0)
	y.Scal(0.0)

	s := matrix.FloatZeros(cdim, 1)
	z := matrix.FloatZeros(cdim, 1)
	dx := c.Copy()
	dy := b.Copy()
	ds := matrix.FloatZeros(cdim, 1)
	dz := matrix.FloatZeros(cdim, 1)
	// these are singleton matrix
	dkappa := matrix.FloatValue(0.0)
	dtau := matrix.FloatValue(0.0)

	checkpnt.AddVerifiable("x", x)
	checkpnt.AddMatrixVar("s", s)
	checkpnt.AddMatrixVar("z", z)
	checkpnt.AddVerifiable("dx", dx)
	checkpnt.AddMatrixVar("ds", ds)
	checkpnt.AddMatrixVar("dz", dz)
	checkpnt.Check("00init", 1)

	var W *sets.FloatMatrixSet
	var f KKTFuncVar
	if primalstart == nil || dualstart == nil {
		// Factor
		//
		//     [ 0   A'  G' ]
		//     [ A   0   0  ].
		//     [ G   0  -I  ]
		//
		W = sets.NewFloatSet("d", "di", "v", "beta", "r", "rti")
		dd := dims.At("l")[0]
		mat := matrix.FloatOnes(dd, 1)
		W.Set("d", mat)
		mat = matrix.FloatOnes(dd, 1)
		W.Set("di", mat)
		dq := len(dims.At("q"))
		W.Set("beta", matrix.FloatOnes(dq, 1))

		for _, n := range dims.At("q") {
			vm := matrix.FloatZeros(n, 1)
			vm.SetIndex(0, 1.0)
			W.Append("v", vm)
		}
		for _, n := range dims.At("s") {
			W.Append("r", matrix.FloatIdentity(n))
			W.Append("rti", matrix.FloatIdentity(n))
		}
		f, err = kktsolver(W)
		if err != nil {
			fmt.Printf("kktsolver error: %s\n", err)
			return
		}
		checkpnt.AddScaleVar(W)
	}

	checkpnt.Check("05init", 5)
	if primalstart == nil {
		// minimize    || G * x - h ||^2
		// subject to  A * x = b
		//
		// by solving
		//
		//     [ 0   A'  G' ]   [ x  ]   [ 0 ]
		//     [ A   0   0  ] * [ dy ] = [ b ].
		//     [ G   0  -I  ]   [ -s ]   [ h ]
		//blas.ScalFloat(x, 0.0)
		//blas.CopyFloat(b, dy)
		checkpnt.MinorPush(5)
		x.Scal(0.0)
		mCopy(b, dy)
		blas.CopyFloat(h, s)

		err = f(x, dy, s)
		if err != nil {
			fmt.Printf("f(x,dy,s): %s\n", err)
			return
		}
		blas.ScalFloat(s, -1.0)
		//fmt.Printf("initial s=\n%v\n", s.ToString("%.5f"))
		checkpnt.MinorPop()
	} else {
		mCopy(&matrixVar{primalstart.At("x")[0]}, x)
		blas.Copy(primalstart.At("s")[0], s)
	}

	// ts = min{ t | s + t*e >= 0 }
	ts, _ := maxStep(s, dims, 0, nil)
	if ts >= 0 && primalstart != nil {
		err = errors.New("initial s is not positive")
		return
	}
	//fmt.Printf("initial ts=%.5f\n", ts)
	checkpnt.Check("10init", 10)

	if dualstart == nil {
		// minimize   || z ||^2
		// subject to G'*z + A'*y + c = 0
		//
		// by solving
		//
		//     [ 0   A'  G' ] [ dx ]   [ -c ]
		//     [ A   0   0  ] [ y  ] = [  0 ].
		//     [ G   0  -I  ] [ z  ]   [  0 ]
		//blas.Copy(c, dx)
		//blas.ScalFloat(dx, -1.0)
		//blas.ScalFloat(y, 0.0)
		checkpnt.MinorPush(10)
		mCopy(c, dx)
		dx.Scal(-1.0)
		y.Scal(0.0)
		blas.ScalFloat(z, 0.0)
		err = f(dx, y, z)
		if err != nil {
			fmt.Printf("f(dx,y,z): %s\n", err)
			return
		}
		//fmt.Printf("initial z=\n%v\n", z.ToString("%.5f"))
		checkpnt.MinorPop()
	} else {
		if len(dualstart.At("y")) > 0 {
			mCopy(&matrixVar{dualstart.At("y")[0]}, y)
		}
		blas.Copy(dualstart.At("z")[0], z)
	}

	// ts = min{ t | z + t*e >= 0 }
	tz, _ := maxStep(z, dims, 0, nil)
	if tz >= 0 && dualstart != nil {
		err = errors.New("initial z is not positive")
		return
	}
	//fmt.Printf("initial tz=%.5f\n", tz)

	nrms := snrm2(s, dims, 0)
	nrmz := snrm2(z, dims, 0)

	gap := 0.0
	pcost := 0.0
	dcost := 0.0
	relgap := 0.0

	checkpnt.Check("20init", 0)

	if primalstart == nil && dualstart == nil {
		gap = sdot(s, z, dims, 0)
		pcost = c.Dot(x)
		dcost = -b.Dot(y) - sdot(h, z, dims, 0)
		if pcost < 0.0 {
			relgap = gap / -pcost
		} else if dcost > 0.0 {
			relgap = gap / dcost
		} else {
			relgap = math.NaN()
		}
		if ts <= 0 && tz < 0 &&
			(gap <= absTolerance || (!math.IsNaN(relgap) && relgap <= relTolerance)) {
			// Constructed initial points happen to be feasible and optimal

			ind := dims.At("l")[0] + dims.Sum("q")
			for _, m := range dims.At("s") {
				symm(s, m, ind)
				symm(z, m, ind)
				ind += m * m
			}

			// rx = A'*y + G'*z + c
			rx := c.Copy()
			Af(y, rx, 1.0, 1.0, la.OptTrans)
			Gf(&matrixVar{z}, rx, 1.0, 1.0, la.OptTrans)
			resx := math.Sqrt(rx.Dot(rx))
			// ry = b - A*x
			ry := b.Copy()
			Af(x, ry, -1.0, -1.0, la.OptNoTrans)
			resy := math.Sqrt(ry.Dot(ry))
			// rz = s + G*x - h
			rz := matrix.FloatZeros(cdim, 1)

			Gf(x, &matrixVar{rz}, 1.0, 0.0, la.OptNoTrans)
			blas.AxpyFloat(s, rz, 1.0)
			blas.AxpyFloat(h, rz, -1.0)
			resz := snrm2(rz, dims, 0)

			pres := math.Max(resy/resy0, resz/resz0)
			dres := resx / resx0
			cx := c.Dot(x)
			by := b.Dot(y)
			hz := sdot(h, z, dims, 0)

			//sol.X = x; sol.Y = y; sol.S = s; sol.Z = z
			sol.Result = sets.NewFloatSet("x", "y", "s", "x")
			sol.Result.Append("x", x.Matrix())
			sol.Result.Append("y", y.Matrix())
			sol.Result.Append("s", s)
			sol.Result.Append("z", z)
			sol.Status = Optimal
			sol.Gap = gap
			sol.RelativeGap = relgap
			sol.PrimalObjective = cx
			sol.DualObjective = -(by + hz)
			sol.PrimalInfeasibility = pres
			sol.DualInfeasibility = dres
			sol.PrimalSlack = -ts
			sol.DualSlack = -tz

			return
		}

		if ts >= -1e-8*math.Max(nrms, 1.0) {
			a := 1.0 + ts
			is := make([]int, 0)
			// indexes s[:dims['l']]
			if dims.At("l")[0] > 0 {
				is = append(is, matrix.MakeIndexSet(0, dims.At("l")[0], 1)...)
			}
			// indexes s[indq[:-1]]
			if len(indq) > 1 {
				is = append(is, indq[:len(indq)-1]...)
			}
			// indexes s[ind:ind+m*m:m+1] (diagonal)
			ind := dims.Sum("l", "q")
			for _, m := range dims.At("s") {
				is = append(is, matrix.MakeIndexSet(ind, ind+m*m, m+1)...)
				ind += m * m
			}
			for _, k := range is {
				s.SetIndex(k, a+s.GetIndex(k))
			}
		}

		if tz >= -1e-8*math.Max(nrmz, 1.0) {
			a := 1.0 + tz
			is := make([]int, 0)
			// indexes z[:dims['l']]
			if dims.At("l")[0] > 0 {
				is = append(is, matrix.MakeIndexSet(0, dims.At("l")[0], 1)...)
			}
			// indexes z[indq[:-1]]
			if len(indq) > 1 {
				is = append(is, indq[:len(indq)-1]...)
			}
			// indexes z[ind:ind+m*m:m+1] (diagonal)
			ind := dims.Sum("l", "q")
			for _, m := range dims.At("s") {
				is = append(is, matrix.MakeIndexSet(ind, ind+m*m, m+1)...)
				ind += m * m
			}
			for _, k := range is {
				z.SetIndex(k, a+z.GetIndex(k))
			}
		}
	} else if primalstart == nil && dualstart != nil {
		if ts >= -1e-8*math.Max(nrms, 1.0) {
			a := 1.0 + ts
			is := make([]int, 0)
			if dims.At("l")[0] > 0 {
				is = append(is, matrix.MakeIndexSet(0, dims.At("l")[0], 1)...)
			}
			if len(indq) > 1 {
				is = append(is, indq[:len(indq)-1]...)
			}
			ind := dims.Sum("l", "q")
			for _, m := range dims.At("s") {
				is = append(is, matrix.MakeIndexSet(ind, ind+m*m, m+1)...)
				ind += m * m
			}
			for _, k := range is {
				s.SetIndex(k, a+s.GetIndex(k))
			}
		}
	} else if primalstart != nil && dualstart == nil {
		if tz >= -1e-8*math.Max(nrmz, 1.0) {
			a := 1.0 + tz
			is := make([]int, 0)
			if dims.At("l")[0] > 0 {
				is = append(is, matrix.MakeIndexSet(0, dims.At("l")[0], 1)...)
			}
			if len(indq) > 1 {
				is = append(is, indq[:len(indq)-1]...)
			}
			ind := dims.Sum("l", "q")
			for _, m := range dims.At("s") {
				is = append(is, matrix.MakeIndexSet(ind, ind+m*m, m+1)...)
				ind += m * m
			}
			for _, k := range is {
				z.SetIndex(k, a+z.GetIndex(k))
			}
		}
	}

	tau := matrix.FloatValue(1.0)
	kappa := matrix.FloatValue(1.0)
	wkappa3 := matrix.FloatValue(0.0)

	rx := c.Copy()
	hrx := c.Copy()
	ry := b.Copy()
	hry := b.Copy()
	rz := matrix.FloatZeros(cdim, 1)
	hrz := matrix.FloatZeros(cdim, 1)
	sigs := matrix.FloatZeros(dims.Sum("s"), 1)
	sigz := matrix.FloatZeros(dims.Sum("s"), 1)
	lmbda := matrix.FloatZeros(cdim_diag+1, 1)
	lmbdasq := matrix.FloatZeros(cdim_diag+1, 1)

	gap = sdot(s, z, dims, 0)

	var x1, y1 MatrixVariable
	var z1 *matrix.FloatMatrix
	var dg, dgi float64
	var th *matrix.FloatMatrix
	var WS fVarClosure
	var f3 KKTFuncVar
	var cx, by, hz, rt float64
	var hresx, resx, hresy, resy, hresz, resz float64
	var dres, pres, dinfres, pinfres float64

	// check point variables
	checkpnt.AddMatrixVar("lmbda", lmbda)
	checkpnt.AddMatrixVar("lmbdasq", lmbdasq)
	checkpnt.AddVerifiable("rx", rx)
	checkpnt.AddVerifiable("ry", ry)
	checkpnt.AddMatrixVar("rz", rz)
	checkpnt.AddFloatVar("resx", &resx)
	checkpnt.AddFloatVar("resy", &resy)
	checkpnt.AddFloatVar("resz", &resz)
	checkpnt.AddFloatVar("hresx", &hresx)
	checkpnt.AddFloatVar("hresy", &hresy)
	checkpnt.AddFloatVar("hresz", &hresz)
	checkpnt.AddFloatVar("cx", &cx)
	checkpnt.AddFloatVar("by", &by)
	checkpnt.AddFloatVar("hz", &hz)
	checkpnt.AddFloatVar("gap", &gap)
	checkpnt.AddFloatVar("pres", &pres)
	checkpnt.AddFloatVar("dres", &dres)

	for iter := 0; iter < maxIter+1; iter++ {
		checkpnt.MajorNext()
		checkpnt.Check("loop-start", 100)

		// hrx = -A'*y - G'*z
		Af(y, hrx, -1.0, 0.0, la.OptTrans)
		Gf(&matrixVar{z}, hrx, -1.0, 1.0, la.OptTrans)
		hresx = math.Sqrt(hrx.Dot(hrx))

		// rx = hrx - c*tau
		//    = -A'*y - G'*z - c*tau
		mCopy(hrx, rx)
		c.Axpy(rx, -tau.Float())
		resx = math.Sqrt(rx.Dot(rx)) / tau.Float()

		// hry = A*x
		Af(x, hry, 1.0, 0.0, la.OptNoTrans)
		hresy = math.Sqrt(hry.Dot(hry))

		// ry = hry - b*tau
		//    = A*x - b*tau
		mCopy(hry, ry)
		b.Axpy(ry, -tau.Float())
		resy = math.Sqrt(ry.Dot(ry)) / tau.Float()

		// hrz = s + G*x
		Gf(x, &matrixVar{hrz}, 1.0, 0.0, la.OptNoTrans)
		blas.AxpyFloat(s, hrz, 1.0)
		hresz = snrm2(hrz, dims, 0)

		// rz = hrz - h*tau
		//    = s + G*x - h*tau
		blas.ScalFloat(rz, 0.0)
		blas.AxpyFloat(hrz, rz, 1.0)
		blas.AxpyFloat(h, rz, -tau.Float())
		resz = snrm2(rz, dims, 0) / tau.Float()

		// rt = kappa + c'*x + b'*y + h'*z '
		cx = c.Dot(x)
		by = b.Dot(y)
		hz = sdot(h, z, dims, 0)
		rt = kappa.Float() + cx + by + hz

		// Statistics for stopping criteria
		pcost = cx / tau.Float()
		dcost = -(by + hz) / tau.Float()

		if pcost < 0.0 {
			relgap = gap / -pcost
		} else if dcost > 0.0 {
			relgap = gap / dcost
		} else {
			relgap = math.NaN()
		}

		pres = math.Max(resy/resy0, resz/resz0)
		dres = resx / resx0
		pinfres = math.NaN()
		if hz+by < 0.0 {
			pinfres = hresx / resx0 / (-hz - by)
		}
		dinfres = math.NaN()
		if cx < 0.0 {
			dinfres = math.Max(hresy/resy0, hresz/resz0) / (-cx)
		}

		if solopts.ShowProgress {
			if iter == 0 {
				// show headers of something
				fmt.Printf("% 10s% 12s% 10s% 8s% 7s % 5s\n",
					"pcost", "dcost", "gap", "pres", "dres", "k/t")
			}
			// show something
			fmt.Printf("%2d: % 8.4e % 8.4e % 4.0e% 7.0e% 7.0e% 7.0e\n",
				iter, pcost, dcost, gap, pres, dres, kappa.GetIndex(0)/tau.GetIndex(0))
		}

		checkpnt.Check("isready", 200)

		if (pres <= feasTolerance && dres <= feasTolerance &&
			(gap <= absTolerance || (!math.IsNaN(relgap) && relgap <= relTolerance))) ||
			iter == maxIter {
			// done
			x.Scal(1.0 / tau.Float())
			y.Scal(1.0 / tau.Float())
			blas.ScalFloat(s, 1.0/tau.Float())
			blas.ScalFloat(z, 1.0/tau.Float())
			ind := dims.Sum("l", "q")
			for _, m := range dims.At("s") {
				symm(s, m, ind)
				symm(z, m, ind)
				ind += m * m
			}
			ts, _ = maxStep(s, dims, 0, nil)
			tz, _ = maxStep(z, dims, 0, nil)
			if iter == maxIter {
				// MaxIterations exceeded
				if solopts.ShowProgress {
					fmt.Printf("No solution. Max iterations exceeded\n")
				}
				err = errors.New("No solution. Max iterations exceeded")
				//sol.X = x; sol.Y = y; sol.S = s; sol.Z = z
				sol.Result = sets.NewFloatSet("x", "y", "s", "x")
				sol.Result.Append("x", x.Matrix())
				sol.Result.Append("y", y.Matrix())
				sol.Result.Append("s", s)
				sol.Result.Append("z", z)
				sol.Status = Unknown
				sol.Gap = gap
				sol.RelativeGap = relgap
				sol.PrimalObjective = pcost
				sol.DualObjective = dcost
				sol.PrimalInfeasibility = pres
				sol.DualInfeasibility = dres
				sol.PrimalSlack = -ts
				sol.DualSlack = -tz
				sol.PrimalResidualCert = pinfres
				sol.DualResidualCert = dinfres
				sol.Iterations = iter
				return
			} else {
				// Optimal
				if solopts.ShowProgress {
					fmt.Printf("Optimal solution.\n")
				}
				err = nil
				//sol.X = x; sol.Y = y; sol.S = s; sol.Z = z
				sol.Result = sets.NewFloatSet("x", "y", "s", "x")
				sol.Result.Append("x", x.Matrix())
				sol.Result.Append("y", y.Matrix())
				sol.Result.Append("s", s)
				sol.Result.Append("z", z)
				sol.Status = Optimal
				sol.Gap = gap
				sol.RelativeGap = relgap
				sol.PrimalObjective = pcost
				sol.DualObjective = dcost
				sol.PrimalInfeasibility = pres
				sol.DualInfeasibility = dres
				sol.PrimalSlack = -ts
				sol.DualSlack = -tz
				sol.PrimalResidualCert = math.NaN()
				sol.DualResidualCert = math.NaN()
				sol.Iterations = iter
				return
			}
		} else if !math.IsNaN(pinfres) && pinfres <= feasTolerance {
			// Primal Infeasible
			if solopts.ShowProgress {
				fmt.Printf("Primal infeasible.\n")
			}
			err = errors.New("Primal infeasible")
			y.Scal(1.0 / (-hz - by))
			blas.ScalFloat(z, 1.0/(-hz-by))
			//sol.X = nil; sol.Y = nil; sol.S = nil; sol.Z = nil
			ind := dims.Sum("l", "q")
			for _, m := range dims.At("s") {
				symm(z, m, ind)
				ind += m * m
			}
			tz, _ = maxStep(z, dims, 0, nil)
			sol.Status = PrimalInfeasible
			sol.Result = sets.NewFloatSet("x", "y", "s", "x")
			sol.Result.Append("x", nil)
			sol.Result.Append("y", nil)
			sol.Result.Append("s", nil)
			sol.Result.Append("z", nil)
			sol.Gap = math.NaN()
			sol.RelativeGap = math.NaN()
			sol.PrimalObjective = math.NaN()
			sol.DualObjective = 1.0
			sol.PrimalInfeasibility = math.NaN()
			sol.DualInfeasibility = math.NaN()
			sol.PrimalSlack = math.NaN()
			sol.DualSlack = -tz
			sol.PrimalResidualCert = pinfres
			sol.DualResidualCert = math.NaN()
			sol.Iterations = iter
			return
		} else if !math.IsNaN(dinfres) && dinfres <= feasTolerance {
			// Dual Infeasible
			if solopts.ShowProgress {
				fmt.Printf("Dual infeasible.\n")
			}
			err = errors.New("Primal infeasible")
			x.Scal(1.0 / (-cx))
			blas.ScalFloat(s, 1.0/(-cx))
			//sol.X = nil; sol.Y = nil; sol.S = nil; sol.Z = nil
			ind := dims.Sum("l", "q")
			for _, m := range dims.At("s") {
				symm(s, m, ind)
				ind += m * m
			}
			ts, _ = maxStep(s, dims, 0, nil)
			sol.Status = PrimalInfeasible
			sol.Result = sets.NewFloatSet("x", "y", "s", "x")
			sol.Result.Append("x", nil)
			sol.Result.Append("y", nil)
			sol.Result.Append("s", nil)
			sol.Result.Append("z", nil)
			sol.Gap = math.NaN()
			sol.RelativeGap = math.NaN()
			sol.PrimalObjective = 1.0
			sol.DualObjective = math.NaN()
			sol.PrimalInfeasibility = math.NaN()
			sol.DualInfeasibility = math.NaN()
			sol.PrimalSlack = -ts
			sol.DualSlack = math.NaN()
			sol.PrimalResidualCert = math.NaN()
			sol.DualResidualCert = dinfres
			sol.Iterations = iter
			return
		}

		// Compute initial scaling W:
		//
		//     W * z = W^{-T} * s = lambda
		//     dg * tau = 1/dg * kappa = lambdag.
		if iter == 0 {
			//fmt.Printf("compute scaling: lmbda=\n%v\n", lmbda.ToString("%.17f"))
			//fmt.Printf("s=\n%v\n", s.ToString("%.17f"))
			//fmt.Printf("z=\n%v\n", z.ToString("%.17f"))
			W, err = computeScaling(s, z, lmbda, dims, 0)
			checkpnt.AddScaleVar(W)

			//     dg = sqrt( kappa / tau )
			//     dgi = sqrt( tau / kappa )
			//     lambda_g = sqrt( tau * kappa )
			//
			// lambda_g is stored in the last position of lmbda.

			dg = math.Sqrt(kappa.Float() / tau.Float())
			dgi = math.Sqrt(float64(tau.Float() / kappa.Float()))
			lmbda.SetIndex(-1, math.Sqrt(float64(tau.Float()*kappa.Float())))
			//fmt.Printf("lmbda=\n%v\n", lmbda.ToString("%.17f"))
			//W.Print()
			checkpnt.Check("compute_scaling", 300)
		}
		// lmbdasq := lmbda o lmbda
		ssqr(lmbdasq, lmbda, dims, 0)
		lmbdasq.SetIndex(-1, lmbda.GetIndex(-1)*lmbda.GetIndex(-1))

		// f3(x, y, z) solves
		//
		//     [ 0  A'  G'   ] [ ux        ]   [ bx ]
		//     [ A  0   0    ] [ uy        ] = [ by ].
		//     [ G  0  -W'*W ] [ W^{-1}*uz ]   [ bz ]
		//
		// On entry, x, y, z contain bx, by, bz.
		// On exit, they contain ux, uy, uz.
		//
		// Also solve
		//
		//     [ 0   A'  G'    ] [ x1        ]          [ c ]
		//     [-A   0   0     ]*[ y1        ] = -dgi * [ b ].
		//     [-G   0   W'*W  ] [ W^{-1}*z1 ]          [ h ]

		f3, err = kktsolver(W)
		if err != nil {
			fmt.Printf("kktsolver error=%v\n", err)
			return
		}
		if iter == 0 {
			x1 = c.Copy()
			y1 = b.Copy()
			z1 = matrix.FloatZeros(cdim, 1)
			checkpnt.AddVerifiable("x1", x1)
			checkpnt.AddMatrixVar("z1", z1)
		}
		mCopy(c, x1)
		x1.Scal(-1.0)
		mCopy(b, y1)
		blas.Copy(h, z1)
		err = f3(x1, y1, z1)
		//fmt.Printf("f3 result: x1=\n%v\nf3 result: z1=\n%v\n", x1, z1)
		x1.Scal(dgi)
		y1.Scal(dgi)
		blas.ScalFloat(z1, dgi)

		if err != nil {
			if iter == 0 && primalstart != nil && dualstart != nil {
				err = errors.New("Rank(A) < p or Rank([G; A]) < n")
				return
			} else {
				t_ := 1.0 / tau.Float()
				x.Scal(t_)
				y.Scal(t_)
				blas.ScalFloat(s, t_)
				blas.ScalFloat(z, t_)
				ind := dims.Sum("l", "q")
				for _, m := range dims.At("s") {
					symm(s, m, ind)
					symm(z, m, ind)
					ind += m * m
				}
				ts, _ = maxStep(s, dims, 0, nil)
				tz, _ = maxStep(z, dims, 0, nil)
				err = errors.New("Terminated (singular KKT matrix).")
				//sol.X = x; sol.Y = y; sol.S = s; sol.Z = z
				sol.Result = sets.NewFloatSet("x", "y", "s", "x")
				sol.Result.Append("x", x.Matrix())
				sol.Result.Append("y", y.Matrix())
				sol.Result.Append("s", s)
				sol.Result.Append("z", z)
				sol.Status = Unknown
				sol.RelativeGap = relgap
				sol.PrimalObjective = pcost
				sol.DualObjective = dcost
				sol.PrimalInfeasibility = pres
				sol.DualInfeasibility = dres
				sol.PrimalSlack = -ts
				sol.DualSlack = -tz
				sol.Iterations = iter
				return
			}
		}

		// f6_no_ir(x, y, z, tau, s, kappa) solves
		//
		//    [ 0         ]   [  0   A'  G'  c ] [ ux        ]    [ bx   ]
		//    [ 0         ]   [ -A   0   0   b ] [ uy        ]    [ by   ]
		//    [ W'*us     ] - [ -G   0   0   h ] [ W^{-1}*uz ] = -[ bz   ]
		//    [ dg*ukappa ]   [ -c' -b' -h'  0 ] [ utau/dg   ]    [ btau ]
		//
		//    lmbda o (uz + us) = -bs
		//    lmbdag * (utau + ukappa) = -bkappa.
		//
		// On entry, x, y, z, tau, s, kappa contain bx, by, bz, btau,
		// bkappa.  On exit, they contain ux, uy, uz, utau, ukappa.

		// th = W^{-T} * h
		if iter == 0 {
			th = matrix.FloatZeros(cdim, 1)
			checkpnt.AddMatrixVar("th", th)
		}

		blas.Copy(h, th)
		scale(th, W, true, true)

		f6_no_ir := func(x, y MatrixVariable, z, tau, s, kappa *matrix.FloatMatrix) (err error) {
			// Solve
			//
			// [  0   A'  G'    0   ] [ ux        ]
			// [ -A   0   0     b   ] [ uy        ]
			// [ -G   0   W'*W  h   ] [ W^{-1}*uz ]
			// [ -c' -b' -h'    k/t ] [ utau/dg   ]
			//
			//   [ bx                    ]
			//   [ by                    ]
			// = [ bz - W'*(lmbda o\ bs) ]
			//   [ btau - bkappa/tau     ]
			//
			// us = -lmbda o\ bs - uz
			// ukappa = -bkappa/lmbdag - utau.

			// First solve
			//
			// [ 0  A' G'   ] [ ux        ]   [  bx                    ]
			// [ A  0  0    ] [ uy        ] = [ -by                    ]
			// [ G  0 -W'*W ] [ W^{-1}*uz ]   [ -bz + W'*(lmbda o\ bs) ]

			minor := checkpnt.MinorTop()
			err = nil
			// y := -y = -by
			y.Scal(-1.0)

			// s := -lmbda o\ s = -lmbda o\ bs
			err = sinv(s, lmbda, dims, 0)
			blas.ScalFloat(s, -1.0)

			// z := -(z + W'*s) = -bz + W'*(lambda o\ bs)
			blas.Copy(s, ws3)
			checkpnt.Check("prescale", minor+5)
			checkpnt.MinorPush(minor + 5)
			err = scale(ws3, W, true, false)
			checkpnt.MinorPop()
			if err != nil {
				fmt.Printf("scale error: %s\n", err)
			}
			blas.AxpyFloat(ws3, z, 1.0)
			blas.ScalFloat(z, -1.0)

			checkpnt.Check("f3-call", minor+20)
			checkpnt.MinorPush(minor + 20)
			err = f3(x, y, z)
			checkpnt.MinorPop()
			checkpnt.Check("f3-return", minor+40)

			// Combine with solution of
			//
			// [ 0   A'  G'    ] [ x1         ]          [ c ]
			// [-A   0   0     ] [ y1         ] = -dgi * [ b ]
			// [-G   0   W'*W  ] [ W^{-1}*dzl ]          [ h ]
			//
			// to satisfy
			//
			// -c'*x - b'*y - h'*W^{-1}*z + dg*tau = btau - bkappa/tau. '

			// , kappa[0] := -kappa[0] / lmbd[-1] = -bkappa / lmbdag
			kap_ := kappa.Float()
			tau_ := tau.Float()
			kap_ = -kap_ / lmbda.GetIndex(-1)
			// tau[0] = tau[0] + kappa[0] / dgi = btau[0] - bkappa / tau
			tau_ = tau_ + kap_/dgi

			//tau[0] = dgi * ( tau[0] + xdot(c,x) + ydot(b,y) +
			//    misc.sdot(th, z, dims) ) / (1.0 + misc.sdot(z1, z1, dims))
			//tau_ = tau_ + blas.DotFloat(c, x) + blas.DotFloat(b, y) + sdot(th, z, dims, 0)
			tau_ += c.Dot(x)
			tau_ += b.Dot(y)
			tau_ += sdot(th, z, dims, 0)
			tau_ = dgi * tau_ / (1.0 + sdot(z1, z1, dims, 0))
			tau.SetValue(tau_)
			x1.Axpy(x, tau_)
			y1.Axpy(y, tau_)
			blas.AxpyFloat(z1, z, tau_)

			blas.AxpyFloat(z, s, -1.0)
			kap_ = kap_ - tau_
			kappa.SetValue(kap_)
			return
		}

		// f6(x, y, z, tau, s, kappa) solves the same system as f6_no_ir,
		// but applies iterative refinement. Following variables part of f6-closure
		// and ~ 12 is the limit. We wrap them to a structure.

		if iter == 0 {
			if refinement > 0 || solopts.Debug {
				WS.wx = c.Copy()
				WS.wy = b.Copy()
				WS.wz = matrix.FloatZeros(cdim, 1)
				WS.ws = matrix.FloatZeros(cdim, 1)
				WS.wtau = matrix.FloatValue(0.0)
				WS.wkappa = matrix.FloatValue(0.0)
				checkpnt.AddVerifiable("wx", WS.wx)
				checkpnt.AddMatrixVar("wz", WS.wz)
				checkpnt.AddMatrixVar("ws", WS.ws)
			}
			if refinement > 0 {
				WS.wx2 = c.Copy()
				WS.wy2 = b.Copy()
				WS.wz2 = matrix.FloatZeros(cdim, 1)
				WS.ws2 = matrix.FloatZeros(cdim, 1)
				WS.wtau2 = matrix.FloatValue(0.0)
				WS.wkappa2 = matrix.FloatValue(0.0)
				checkpnt.AddVerifiable("wx2", WS.wx2)
				checkpnt.AddMatrixVar("wz2", WS.wz2)
				checkpnt.AddMatrixVar("ws2", WS.ws2)
			}
		}

		f6 := func(x, y MatrixVariable, z, tau, s, kappa *matrix.FloatMatrix) error {
			var err error = nil
			minor := checkpnt.MinorTop()
			checkpnt.Check("startf6", minor+100)
			if refinement > 0 || solopts.Debug {
				mCopy(x, WS.wx)
				mCopy(y, WS.wy)
				blas.Copy(z, WS.wz)
				blas.Copy(s, WS.ws)
				WS.wtau.SetValue(tau.Float())
				WS.wkappa.SetValue(kappa.Float())
			}
			checkpnt.Check("pref6_no_ir", minor+200)
			err = f6_no_ir(x, y, z, tau, s, kappa)
			checkpnt.Check("postf6_no_ir", minor+399)
			for i := 0; i < refinement; i++ {
				mCopy(WS.wx, WS.wx2)
				mCopy(WS.wy, WS.wy2)
				blas.Copy(WS.wz, WS.wz2)
				blas.Copy(WS.ws, WS.ws2)
				WS.wtau2.SetValue(WS.wtau.Float())
				WS.wkappa2.SetValue(WS.wkappa.Float())
				checkpnt.Check("res-call", minor+400)

				checkpnt.MinorPush(minor + 400)
				err = res(x, y, z, tau, s, kappa, WS.wx2, WS.wy2, WS.wz2, WS.wtau2, WS.ws2, WS.wkappa2, W, dg, lmbda)
				checkpnt.MinorPop()

				checkpnt.Check("refine_pref6_no_ir", minor+500)
				checkpnt.MinorPush(minor + 500)
				err = f6_no_ir(WS.wx2, WS.wy2, WS.wz2, WS.wtau2, WS.ws2, WS.wkappa2)
				checkpnt.MinorPop()

				checkpnt.Check("refine_postf6_no_ir", minor+600)
				WS.wx2.Axpy(x, 1.0)
				WS.wy2.Axpy(y, 1.0)
				blas.AxpyFloat(WS.wz2, z, 1.0)
				blas.AxpyFloat(WS.ws2, s, 1.0)
				tau.SetValue(tau.Float() + WS.wtau2.Float())
				kappa.SetValue(kappa.Float() + WS.wkappa2.Float())
			}
			if solopts.Debug {
				checkpnt.MinorPush(minor + 700)
				res(x, y, z, tau, s, kappa, WS.wx, WS.wy, WS.wz, WS.wtau, WS.ws, WS.wkappa, W, dg, lmbda)
				checkpnt.MinorPop()
				fmt.Printf("KKT residuals\n")
				fmt.Printf("    'x'    : %.6e\n", math.Sqrt(WS.wx.Dot(WS.wx)))
				fmt.Printf("    'y'    : %.6e\n", math.Sqrt(WS.wy.Dot(WS.wy)))
				fmt.Printf("    'z'    : %.6e\n", snrm2(WS.wz, dims, 0))
				fmt.Printf("    'tau'  : %.6e\n", math.Abs(WS.wtau.Float()))
				fmt.Printf("    's'    : %.6e\n", snrm2(WS.ws, dims, 0))
				fmt.Printf("    'kappa': %.6e\n", math.Abs(WS.wkappa.Float()))
			}
			return err
		}

		var nrm float64 = blas.Nrm2Float(lmbda)
		mu := math.Pow(nrm, 2.0) / (1.0 + float64(cdim_diag))
		sigma := 0.0
		var step, tt, tk float64

		for i := 0; i < 2; i++ {
			// Solve
			//
			// [ 0         ]   [  0   A'  G'  c ] [ dx        ]
			// [ 0         ]   [ -A   0   0   b ] [ dy        ]
			// [ W'*ds     ] - [ -G   0   0   h ] [ W^{-1}*dz ]
			// [ dg*dkappa ]   [ -c' -b' -h'  0 ] [ dtau/dg   ]
			//
			//               [ rx   ]
			//               [ ry   ]
			// = - (1-sigma) [ rz   ]
			//               [ rtau ]
			//
			// lmbda o (dz + ds) = -lmbda o lmbda + sigma*mu*e
			// lmbdag * (dtau + dkappa) = - kappa * tau + sigma*mu
			//
			// ds = -lmbdasq if i is 0
			//    = -lmbdasq - dsa o dza + sigma*mu*e if i is 1
			// dkappa = -lambdasq[-1] if i is 0
			//        = -lambdasq[-1] - dkappaa*dtaua + sigma*mu if i is 1.
			ind := dims.Sum("l", "q")
			ind2 := ind
			blas.Copy(lmbdasq, ds, &la.IOpt{"n", ind})
			blas.ScalFloat(ds, 0.0, &la.IOpt{"offset", ind})
			for _, m := range dims.At("s") {
				blas.Copy(lmbdasq, ds, &la.IOpt{"n", m}, &la.IOpt{"offsetx", ind2},
					&la.IOpt{"offsety", ind}, &la.IOpt{"incy", m + 1})
				ind += m * m
				ind2 += m
			}
			// dkappa[0] = lmbdasq[-1]
			dkappa.SetValue(lmbdasq.GetIndex(-1))

			if i == 1 {
				blas.AxpyFloat(ws3, ds, 1.0)
				ind = dims.Sum("l", "q")
				is := make([]int, 0)
				// indexes: [:dims['l']]
				if dims.At("l")[0] > 0 {
					is = append(is, matrix.MakeIndexSet(0, dims.At("l")[0], 1)...)
				}
				// ...[indq[:-1]]
				if len(indq) > 1 {
					is = append(is, indq[:len(indq)-1]...)
				}
				// ...[ind : ind+m*m : m+1] (diagonal)
				for _, m := range dims.At("s") {
					is = append(is, matrix.MakeIndexSet(ind, ind+m*m, m+1)...)
					ind += m * m
				}
				//ds.Add(-sigma*mu, is...)
				for _, k := range is {
					ds.SetIndex(k, ds.GetIndex(k)-sigma*mu)
				}

				dk_ := dkappa.Float()
				wk_ := wkappa3.Float()
				dkappa.SetValue(dk_ + wk_ - sigma*mu)
			}
			// (dx, dy, dz, dtau) = (1-sigma)*(rx, ry, rz, rt)
			mCopy(rx, dx)
			dx.Scal(1.0 - sigma)
			mCopy(ry, dy)
			dy.Scal(1.0 - sigma)
			blas.Copy(rz, dz)
			blas.ScalFloat(dz, 1.0-sigma)
			// dtau[0] = (1.0 - sigma) * rt
			dtau.SetValue((1.0 - sigma) * rt)

			checkpnt.Check("pref6", (1+i)*1000)
			checkpnt.MinorPush((1 + i) * 1000)
			err = f6(dx, dy, dz, dtau, ds, dkappa)
			checkpnt.MinorPop()
			checkpnt.Check("postf6", (1+i)*1000+800)

			// Save ds o dz and dkappa * dtau for Mehrotra correction
			if i == 0 {
				blas.Copy(ds, ws3)
				sprod(ws3, dz, dims, 0)
				wkappa3.SetValue(dtau.Float() * dkappa.Float())
			}

			// Maximum step to boundary.
			//
			// If i is 1, also compute eigenvalue decomposition of the 's'
			// blocks in ds, dz.  The eigenvectors Qs, Qz are stored in
			// dsk, dzk.  The eigenvalues are stored in sigs, sigz.
			var ts, tz float64

			checkpnt.MinorPush((1+i)*1000 + 900)
			scale2(lmbda, ds, dims, 0, false)
			scale2(lmbda, dz, dims, 0, false)
			checkpnt.MinorPop()
			checkpnt.Check("post-scale2", (1+i)*1000+990)
			if i == 0 {
				ts, _ = maxStep(ds, dims, 0, nil)
				tz, _ = maxStep(dz, dims, 0, nil)
			} else {
				ts, _ = maxStep(ds, dims, 0, sigs)
				tz, _ = maxStep(dz, dims, 0, sigz)
			}
			dt_ := dtau.Float()
			dk_ := dkappa.Float()
			tt = -dt_ / lmbda.GetIndex(-1)
			tk = -dk_ / lmbda.GetIndex(-1)
			t := maxvec([]float64{0.0, ts, tz, tt, tk})
			if t == 0.0 {
				step = 1.0
			} else {
				if i == 0 {
					step = math.Min(1.0, 1.0/t)
				} else {
					step = math.Min(1.0, STEP/t)
				}
			}
			if i == 0 {
				// sigma = (1 - step)^3
				sigma = (1.0 - step) * (1.0 - step) * (1.0 - step)
				//sigma = math.Pow((1.0 - step), EXPON)
			}
		}
		//fmt.Printf("** tau = %.17f, kappa = %.17f\n", tau.Float(), kappa.Float())
		//fmt.Printf("** step = %.17f, sigma = %.17f\n", step, sigma)

		checkpnt.Check("update-xy", 7000)
		// Update x, y
		dx.Axpy(x, step)
		dy.Axpy(y, step)

		// Replace 'l' and 'q' blocks of ds and dz with the updated
		// variables in the current scaling.
		// Replace 's' blocks of ds and dz with the factors Ls, Lz in a
		// factorization Ls*Ls', Lz*Lz' of the updated variables in the
		// current scaling.
		//
		// ds := e + step*ds for 'l' and 'q' blocks.
		// dz := e + step*dz for 'l' and 'q' blocks.
		blas.ScalFloat(ds, step, &la.IOpt{"n", dims.Sum("l", "q")})
		blas.ScalFloat(dz, step, &la.IOpt{"n", dims.Sum("l", "q")})

		is := make([]int, 0)
		is = append(is, matrix.MakeIndexSet(0, dims.At("l")[0], 1)...)
		is = append(is, indq[:len(indq)-1]...)
		for _, k := range is {
			ds.SetIndex(k, 1.0+ds.GetIndex(k))
			dz.SetIndex(k, 1.0+dz.GetIndex(k))
		}
		checkpnt.Check("update-dsdz", 7500)

		// ds := H(lambda)^{-1/2} * ds and dz := H(lambda)^{-1/2} * dz.
		//
		// This replaces the 'l' and 'q' components of ds and dz with the
		// updated variables in the current scaling.
		// The 's' components of ds and dz are replaced with
		//
		// diag(lmbda_k)^{1/2} * Qs * diag(lmbda_k)^{1/2}
		// diag(lmbda_k)^{1/2} * Qz * diag(lmbda_k)^{1/2}
		checkpnt.MinorPush(7500)
		scale2(lmbda, ds, dims, 0, true)
		scale2(lmbda, dz, dims, 0, true)
		checkpnt.MinorPop()

		// sigs := ( e + step*sigs ) ./ lambda for 's' blocks.
		// sigz := ( e + step*sigz ) ./ lambda for 's' blocks.
		blas.ScalFloat(sigs, step)
		blas.ScalFloat(sigz, step)
		sigs.Add(1.0)
		sigz.Add(1.0)
		sdimsum := dims.Sum("s")
		qdimsum := dims.Sum("l", "q")
		blas.TbsvFloat(lmbda, sigs, &la.IOpt{"n", sdimsum}, &la.IOpt{"k", 0},
			&la.IOpt{"lda", 1}, &la.IOpt{"offseta", qdimsum})
		blas.TbsvFloat(lmbda, sigz, &la.IOpt{"n", sdimsum}, &la.IOpt{"k", 0},
			&la.IOpt{"lda", 1}, &la.IOpt{"offseta", qdimsum})

		ind2 := qdimsum
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

		checkpnt.Check("pre-update-scaling", 7700)
		err = updateScaling(W, lmbda, ds, dz)
		checkpnt.Check("post-update-scaling", 7800)

		// For kappa, tau block:
		//
		//     dg := sqrt( (kappa + step*dkappa) / (tau + step*dtau) )
		//         = dg * sqrt( (1 - step*tk) / (1 - step*tt) )
		//
		//     lmbda[-1] := sqrt((tau + step*dtau) * (kappa + step*dkappa))
		//                = lmbda[-1] * sqrt(( 1 - step*tt) * (1 - step*tk))
		dg *= math.Sqrt(1.0-step*tk) / math.Sqrt(1.0-step*tt)
		dgi = 1.0 / dg
		a := math.Sqrt(1.0-step*tk) * math.Sqrt(1.0-step*tt)
		lmbda.SetIndex(-1, a*lmbda.GetIndex(-1))

		// Unscale s, z, tau, kappa (unscaled variables are used only to
		// compute feasibility residuals).
		ind := dims.Sum("l", "q")
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

		ind = dims.Sum("l", "q")
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

		kappa.SetValue(lmbda.GetIndex(-1) / dgi)
		tau.SetValue(lmbda.GetIndex(-1) * dgi)
		g := blas.Nrm2Float(lmbda, &la.IOpt{"n", lmbda.Rows() - 1}) / tau.Float()
		gap = g * g
		checkpnt.Check("end-of-loop", 8000)
		//fmt.Printf(" ** kappa=%.10f, tau=%.10f, gap=%.10f\n", kappa.Float(), tau.Float(), gap)

	}
	return
}

// Local Variables:
// tab-width: 4
// indent-tabs-mode: nil
// End:
