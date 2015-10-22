// Copyright (c) Harri Rautila, 2012

// This file is part of github.com/hrautila/cvx package.
// It is free software, distributed under the terms of GNU Lesser General Public
// License Version 3, or any later version. See the COPYING tile included in this archive.

package cvx

import (
	"errors"
	"fmt"
	"github.com/henrylee2cn/algorithm/cvx/sets"
	la "github.com/henrylee2cn/algorithm/linalg"
	"github.com/henrylee2cn/algorithm/linalg/blas"
	"github.com/henrylee2cn/algorithm/matrix"
	"strconv"
	"strings"
)

// Implements MatrixVariable interface for for CP problem variables.
type epigraph struct {
	mt *matrix.FloatMatrix
	tv float64
}

func newEpigraph(m *matrix.FloatMatrix, t float64) *epigraph {
	return &epigraph{m.Copy(), t}
}

func (u *epigraph) Matrix() *matrix.FloatMatrix {
	return u.m()
}

func (u *epigraph) Dot(vArg MatrixVariable) float64 {
	if v, ok := vArg.(*epigraph); ok {
		return blas.DotFloat(u.m(), v.m()) + u.t()*v.t()
	}
	// really this??
	return blas.DotFloat(u.m(), vArg.Matrix())
}

func (u *epigraph) Scal(alpha float64) {
	blas.ScalFloat(u.m(), alpha)
	u.set(u.t() * alpha)
}

func (u *epigraph) Copy() MatrixVariable {
	return newEpigraph(u.m(), u.t())
}

func (u *epigraph) Axpy(vArg MatrixVariable, alpha float64) {
	if v, ok := vArg.(*epigraph); ok {
		blas.AxpyFloat(u.m(), v.m(), alpha)
		v.set(u.t()*alpha + v.t())
	} else {
		//fmt.Printf("epigraph.Axpy() with non-epigraph argument...\n")
		blas.AxpyFloat(u.m(), v.Matrix(), alpha)
	}
}

func (u *epigraph) String() string {
	s := fmt.Sprintf("%v\n", u.m().ToString("%.7f"))
	s += fmt.Sprintf("/%.7f/", u.t())
	return s
}

// Implement Verifiable interface
func (u *epigraph) Verify(line string) float64 {
	diff := 0.0
	lpar := strings.Index(line, "[")
	rpar := strings.LastIndex(line, "]")
	mstart := strings.Index(line, "{")
	refval, _ := matrix.FloatParse(line[mstart:rpar])
	fval, _ := strconv.ParseFloat(strings.Trim(line[lpar+1:mstart], " "), 64)
	diff += blas.Nrm2Float(matrix.Minus(refval, u.m()))
	d := u.t() - fval
	diff += d * d
	return diff
}

func (u *epigraph) ShowError(line string) {
	lpar := strings.Index(line, "[")
	rpar := strings.LastIndex(line, "]")
	mstart := strings.Index(line, "{")
	refval, _ := matrix.FloatParse(line[mstart:rpar])
	fval, _ := strconv.ParseFloat(strings.Trim(line[lpar+1:mstart], " "), 64)
	df := matrix.Minus(u.m(), refval)
	em, _ := matrix.FloatMatrixStacked(matrix.StackRight, u.m(), refval, df)
	fmt.Printf("my data | ref.data | diff \n%v\n", em.ToString("%.4e"))
	d := u.t() - fval
	fmt.Printf("/%.4e/ | /%.4e/ | /%4e/\n", u.t(), fval, d)
}

func (u *epigraph) m() *matrix.FloatMatrix {
	//return (*u)[0]
	return u.mt
}

func (u *epigraph) t() float64 {
	//return (*u)[1].Float()
	return u.tv
}

func (u *epigraph) set(v float64) {
	//(*u)[1].SetIndex(0, v)
	u.tv = v
}

// internal interface to map ConvexProg for epigraph version for CP solver.
type cpConvexProg interface {
	// Returns (mnl, x0) where mln number of nonlinear inequality constraints
	// and x0 is a point in the domain of f.
	F0() (mnl int, x0 MatrixVariable, err error)

	// Returns a tuple (f, Df) where f is of size (mnl, 1) containing f(x)
	// Df is matrix of size (mnl, n) containing the derivatives of f at x:
	// Df[k,:] is the transpose of the gradient of fk at x. If x is not in
	// domf, return non-nil error.
	F1(x MatrixVariable) (f MatrixVariable, Df MatrixVarDf, err error)

	// F(x, z) with z a positive  matrix of size (mnl, 1). Return a tuple
	// (f, Df, H), where f, Df as above. H is matrix of size (n, n).
	F2(x MatrixVariable, z *matrix.FloatMatrix) (f MatrixVariable, Df MatrixVarDf, H MatrixVarH, err error)
}

// internal interface to map ConvexProg for epigraph version for CP solver.
type cpProg struct {
	convexF ConvexProg
}

func (F *cpProg) F0() (mnl int, x0 MatrixVariable, err error) {
	var m0 *matrix.FloatMatrix
	mnl, m0, err = F.convexF.F0()
	if err != nil {
		return
	}
	mnl += 1
	x0 = newEpigraph(m0, 0.0)
	return
}

func (F *cpProg) F1(xa MatrixVariable) (f MatrixVariable, Df MatrixVarDf, err error) {
	f = nil
	Df = nil
	err = nil
	x, x_ok := xa.(*epigraph)
	if !x_ok {
		err = errors.New("'x' argument not an epigraph")
	}
	var f0, Df0 *matrix.FloatMatrix
	f0, Df0, err = F.convexF.F1(x.m())
	if err != nil {
		fmt.Printf("'cpProg.F1' error: %v\n%s\n", err, x)
		return
	}
	f0.Add(-x.t(), 0)
	f = &matrixVar{f0}
	Df = &epigraphDf{Df0}
	return
}

func (F *cpProg) F2(xa MatrixVariable, z *matrix.FloatMatrix) (f MatrixVariable, Df MatrixVarDf, H MatrixVarH, err error) {
	f = nil
	Df = nil
	H = nil
	err = nil
	x, x_ok := xa.(*epigraph)
	if !x_ok {
		err = errors.New("'x' argument not an epigraph")
		return
	}
	var f0, Df0, H0 *matrix.FloatMatrix
	f0, Df0, H0, err = F.convexF.F2(x.m(), z)
	if err != nil {
		return
	}
	f0.Add(-x.t(), 0)
	f = &matrixVar{f0}
	Df = &epigraphDf{Df0}
	H = &epigraphH{H0}
	return
}

// Implement MatrixVarDf interface for standard matrix valued Df.
type epigraphDf struct {
	df *matrix.FloatMatrix
}

func (d *epigraphDf) Df(u, v MatrixVariable, alpha, beta float64, trans la.Option) (err error) {
	err = nil
	if trans.Equal(la.OptNoTrans) {
		u_e, u_ok := u.(*epigraph)
		v_e := v.Matrix()
		if !u_ok {
			fmt.Printf("Df: not a epigraph\n")
			return errors.New("'u' not a epigraph")
		}
		err = blas.GemvFloat(d.df, u_e.m(), v_e, alpha, beta, la.OptNoTrans)
		v_e.Add(-alpha*u_e.t(), 0)
	} else {
		v_e, v_ok := v.(*epigraph)
		u_e := u.Matrix()
		if !v_ok {
			fmt.Printf("Df: not a epigraph\n")
			return errors.New("'v' not a epigraph")
		}
		err = blas.GemvFloat(d.df, u_e, v_e.m(), alpha, beta, la.OptTrans)
		v_e.set(-alpha*u_e.GetIndex(0) + beta*v_e.t())
	}
	return
}

// Implement MatrixVarH interface for standard matrix valued H in CP problems.
type epigraphH struct {
	h *matrix.FloatMatrix
}

func (g *epigraphH) Hf(u, v MatrixVariable, alpha, beta float64) (err error) {
	err = nil
	u_e, u_ok := u.(*epigraph)
	v_e, v_ok := v.(*epigraph)
	if !u_ok {
		err = errors.New("'u' not a epigraph")
		return
	}
	if !v_ok {
		err = errors.New("'v' not a epigraph")
		return
	}
	err = blas.SymvFloat(g.h, u_e.m(), v_e.m(), alpha, beta)
	v_e.set(v_e.t() + beta*v_e.t())

	return
}

// Implement MatrixVarG interface for CP epigraph parameters
type epMatrixG struct {
	G    *matrix.FloatMatrix
	dims *sets.DimensionSet
}

func (g *epMatrixG) Gf(u, v MatrixVariable, alpha, beta float64, trans la.Option) (err error) {

	err = nil
	if trans.Equal(la.OptNoTrans) {
		ue, u_ok := u.(*epigraph)
		ve := v.Matrix()
		if !u_ok {
			fmt.Printf("Df: not a epigraph\n")
			return errors.New("'u' not a epigraph")
		}
		err = sgemv(g.G, ue.m(), ve, alpha, beta, g.dims, trans)
	} else {
		ve, v_ok := v.(*epigraph)
		ue := u.Matrix()
		if !v_ok {
			return errors.New("'v' not a epigraph")
		}
		err = sgemv(g.G, ue, ve.m(), alpha, beta, g.dims, trans)
		ve.set(beta * ve.t())
	}
	return
}

// Implement MatrixVarG interface for CP epigraph parameters over MatrixG variable
type epiMatrixG struct {
	G    MatrixG
	dims *sets.DimensionSet
}

func (g *epiMatrixG) Gf(u, v MatrixVariable, alpha, beta float64, trans la.Option) (err error) {

	err = nil
	if trans.Equal(la.OptNoTrans) {
		ue, u_ok := u.(*epigraph)
		ve := v.Matrix()
		if !u_ok {
			return errors.New("'u' not a epigraph")
		}
		err = g.G.Gf(ue.m(), ve, alpha, beta, trans)
	} else {
		ve, v_ok := v.(*epigraph)
		ue := u.Matrix()
		if !v_ok {
			return errors.New("'v' not a epigraph")
		}
		err = g.G.Gf(ue, ve.m(), alpha, beta, trans)
		ve.set(beta * ve.t())
	}
	return
}

// Implements MatrixVarA interface for CP epigraph problems with matrix valued A
type epMatrixA struct {
	A *matrix.FloatMatrix
}

func (a *epMatrixA) Af(u, v MatrixVariable, alpha, beta float64, trans la.Option) (err error) {
	err = nil
	if trans.Equal(la.OptNoTrans) {
		ue, u_ok := u.(*epigraph)
		ve := v.Matrix()
		if !u_ok {
			return errors.New("'u' not a epigraph")
		}
		err = blas.GemvFloat(a.A, ue.m(), ve, alpha, beta)
	} else {
		ve, v_ok := v.(*epigraph)
		ue := u.Matrix()
		if !v_ok {
			return errors.New("'v' not a epigraph")
		}
		err = blas.GemvFloat(a.A, ue, ve.m(), alpha, beta, trans)
		ve.set(ve.t() * beta)
	}
	return
}

// Implements MatrixVarA interface for CP epigraph problems with interface MatrixA
type epiMatrixA struct {
	A MatrixA
}

func (a *epiMatrixA) Af(u, v MatrixVariable, alpha, beta float64, trans la.Option) (err error) {
	err = nil
	if trans.Equal(la.OptNoTrans) {
		ue, u_ok := u.(*epigraph)
		ve := v.Matrix()
		if !u_ok {
			return errors.New("'u' not a epigraph")
		}
		err = a.A.Af(ue.m(), ve, alpha, beta, trans)
	} else {
		ve, v_ok := v.(*epigraph)
		ue := u.Matrix()
		if !v_ok {
			return errors.New("'v' not a epigraph")
		}
		err = a.A.Af(ue, ve.m(), alpha, beta, trans)
		ve.set(ve.t() * beta)
	}
	return
}

// Solves a convex optimization problem with a linear objective
//
//       minimize    f0(x)
//       subject to  fk(x) <= 0, k = 1, ..., mnl
//                   G*x   <= h
//                   A*x    = b.
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
func Cp(F ConvexProg, G, h, A, b *matrix.FloatMatrix, dims *sets.DimensionSet, solopts *SolverOptions) (sol *Solution, err error) {

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
	if err = checkConeLpDimensions(dims); err != nil {
		return
	}
	cdim := dims.Sum("l", "q") + dims.SumSquared("s")

	if h.Rows() != cdim {
		err = errors.New(fmt.Sprintf("'h' must be float matrix of size (%d,1)", cdim))
		return
	}

	if G == nil {
		G = matrix.FloatZeros(0, x0.Rows())
	}
	if !G.SizeMatch(cdim, x0.Rows()) {
		estr := fmt.Sprintf("'G' must be of size (%d,%d)", cdim, x0.Rows())
		err = errors.New(estr)
		return
	}

	// Check A and set defaults if it is nil
	if A == nil {
		// zeros rows reduces Gemv to vector products
		A = matrix.FloatZeros(0, x0.Rows())
	}
	if A.Cols() != x0.Rows() {
		estr := fmt.Sprintf("'A' must have %d columns", x0.Rows())
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

	solvername := solopts.KKTSolverName
	if len(solvername) == 0 {
		if len(dims.At("q")) > 0 || len(dims.At("s")) > 0 {
			solvername = "chol"
		} else {
			solvername = "chol2"
		}
	}

	c_e := newEpigraph(x0, 1.0)
	blas.ScalFloat(c_e.m(), 0.0)
	//F_e := &cpProg{F}
	G_e := epMatrixG{G, dims}
	A_e := epMatrixA{A}
	b_e := matrixVar{b}

	var factor kktFactor
	var kktsolver KKTCpSolver = nil
	if kktfunc, ok := solvers[solvername]; ok {
		// kkt function returns us problem spesific factor function.
		factor, err = kktfunc(G, dims, A, mnl)
		if err != nil {
			return nil, err
		}
		// solver is
		kktsolver = func(W *sets.FloatMatrixSet, x, z *matrix.FloatMatrix) (KKTFunc, error) {
			_, Df, H, err := F.F2(x, z)
			if err != nil {
				return nil, err
			}
			return factor(W, H, Df.GetSubMatrix(1, 0))
		}
	} else {
		err = errors.New(fmt.Sprintf("solver '%s' not known", solvername))
		return
	}

	return cp_problem(F, c_e, &G_e, h, &A_e, &b_e, dims, kktsolver, solopts, x0, mnl)
}

// Solves a convex optimization problem with a linear objective
//
//       minimize    f0(x)
//       subject to  fk(x) <= 0, k = 1, ..., mnl
//                   G*x   <= h
//                   A*x    = b.
//
// using custom solver for KKT equations.
//
func CpCustomKKT(F ConvexProg, G, h, A, b *matrix.FloatMatrix, dims *sets.DimensionSet,
	kktsolver KKTCpSolver, solopts *SolverOptions) (sol *Solution, err error) {

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
	if err = checkConeLpDimensions(dims); err != nil {
		return
	}
	cdim := dims.Sum("l", "q") + dims.SumSquared("s")

	if h.Rows() != cdim {
		err = errors.New(fmt.Sprintf("'h' must be float matrix of size (%d,1)", cdim))
		return
	}

	if G == nil {
		G = matrix.FloatZeros(0, x0.Rows())
	}
	if !G.SizeMatch(cdim, x0.Rows()) {
		estr := fmt.Sprintf("'G' must be of size (%d,%d)", cdim, x0.Rows())
		err = errors.New(estr)
		return
	}

	// Check A and set defaults if it is nil
	if A == nil {
		// zeros rows reduces Gemv to vector products
		A = matrix.FloatZeros(0, x0.Rows())
	}
	if A.Cols() != x0.Rows() {
		estr := fmt.Sprintf("'A' must have %d columns", x0.Rows())
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

	if kktsolver == nil {
		err = errors.New("'kktsolver' must be non-nil function.")
		return
	}

	c_e := newEpigraph(x0, 1.0)
	blas.ScalFloat(x0, 0.0)
	G_e := epMatrixG{G, dims}
	A_e := epMatrixA{A}
	b_e := matrixVar{b}

	return cp_problem(F, c_e, &G_e, h, &A_e, &b_e, dims, kktsolver, solopts, x0, mnl)
}

// Solves a convex optimization problem with a linear objective
//
//       minimize    f0(x)
//       subject to  fk(x) <= 0, k = 1, ..., mnl
//                   G*x   <= h
//                   A*x    = b.
//
// using custom solver for KKT equations and constraint equations G and A.
//
func CpCustomMatrix(F ConvexProg, G MatrixG, h *matrix.FloatMatrix, A MatrixA,
	b *matrix.FloatMatrix, dims *sets.DimensionSet, kktsolver KKTCpSolver,
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
	if err = checkConeLpDimensions(dims); err != nil {
		return
	}
	cdim := dims.Sum("l", "q") + dims.SumSquared("s")

	if h.Rows() != cdim {
		err = errors.New(fmt.Sprintf("'h' must be float matrix of size (%d,1)", cdim))
		return
	}

	var G_e MatrixVarG = nil
	if G == nil {
		G_e = &epMatrixG{matrix.FloatZeros(0, x0.Rows()), dims}
	} else {
		G_e = &epiMatrixG{G, dims}
	}

	var A_e MatrixVarA = nil
	if A == nil {
		A_e = &epMatrixA{matrix.FloatZeros(0, x0.Rows())}
	} else {
		A_e = &epiMatrixA{A}
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

	if kktsolver == nil {
		err = errors.New("'kktsolver' must be non-nil function.")
		return
	}

	c_e := newEpigraph(x0, 1.0)
	blas.ScalFloat(c_e.m(), 0.0)
	b_e := matrixVar{b}

	return cp_problem(F, c_e, G_e, h, A_e, &b_e, dims, kktsolver, solopts, x0, mnl)

}

// Here problem is already translated to epigraph format except original convex problem.
// We wrap it and create special CP epigraph kktsolver.
func cp_problem(F ConvexProg, c MatrixVariable, G MatrixVarG, h *matrix.FloatMatrix, A MatrixVarA,
	b MatrixVariable, dims *sets.DimensionSet, kktsolver KKTCpSolver,
	solopts *SolverOptions, x0 *matrix.FloatMatrix, mnl int) (sol *Solution, err error) {

	err = nil

	F_e := &cpProg{F}

	//mx0 := newEpigraph(x0, 0.0)
	cdim := dims.Sum("l", "q") + dims.SumSquared("s")
	ux := x0.Copy()
	uz := matrix.FloatZeros(mnl+cdim, 1)

	kktsolver_e := func(W *sets.FloatMatrixSet, xa MatrixVariable, znl *matrix.FloatMatrix) (KKTFuncVar, error) {
		x, x_ok := xa.(*epigraph)
		_ = x_ok
		We := W.Copy()
		// dnl is matrix
		dnl := W.At("dnl")[0]
		dnli := W.At("dnli")[0]
		We.Set("dnl", matrix.FloatVector(dnl.FloatArray()[1:]))
		We.Set("dnli", matrix.FloatVector(dnli.FloatArray()[1:]))
		g, err := kktsolver(We, x.m(), znl)
		_, Df, _ := F.F1(x.m())
		gradf0 := Df.GetRow(0, nil).Transpose()

		solve := func(xa, ya MatrixVariable, z *matrix.FloatMatrix) (err error) {
			x, x_ok := xa.(*epigraph)
			_ = x_ok // TODO: remove or use x_ok
			y := ya.Matrix()
			err = nil
			a := z.GetIndex(0)
			blas.Copy(x.m(), ux)
			blas.AxpyFloat(gradf0, ux, x.t())
			blas.Copy(z, uz, &la.IOpt{"offsetx", 1})
			err = g(ux, y, uz)
			z.SetIndex(0, -x.t()*dnl.GetIndex(0))
			blas.Copy(uz, z, &la.IOpt{"offsety", 1})
			blas.Copy(ux, x.m())
			val := blas.DotFloat(gradf0, x.m()) + dnl.GetIndex(0)*dnl.GetIndex(0)*x.t() - a
			x.set(val)
			return
		}
		return solve, err
	}
	return cpl_solver(F_e, c, G, h, A, b, dims, kktsolver_e, solopts, nil, mnl)
}

// Local Variables:
// tab-width: 4
// End:
