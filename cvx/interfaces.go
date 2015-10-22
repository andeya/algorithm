// Copyright (c) Harri Rautila, 2012

// This file is part of github.com/hrautila/cvx package.
// It is free software, distributed under the terms of GNU Lesser General Public
// License Version 3, or any later version. See the COPYING tile included in this archive.

package cvx

import (
	"fmt"
	"github.com/henrylee2cn/algorithm/cvx/sets"
	la "github.com/henrylee2cn/algorithm/linalg"
	"github.com/henrylee2cn/algorithm/linalg/blas"
	"github.com/henrylee2cn/algorithm/matrix"
)

// Public interface to provide custom G matrix-vector products
//
// The call Gf(u, v, alpha, beta, trans) should evaluate the matrix-vector products
//
//   v := alpha * G * u + beta * v   if trans is linalg.OptNoTrans
//   v := alpha * G' * u + beta * v  if trans is linalg.OptTrans
//
type MatrixG interface {
	Gf(u, v *matrix.FloatMatrix, alpha, beta float64, trans la.Option) error
}

//
// MatrixVarG provides interface for extended customization with non-matrix type
// primal and dual variables.
//
type MatrixVarG interface {
	Gf(u, v MatrixVariable, alpha, beta float64, trans la.Option) error
}

// Public interface to provide custom A matrix-vector products
//
// The call Af(u, v, alpha, beta, trans) should evaluate the matrix-vector products
//
//   v := alpha * A * u + beta * v   if trans is linalg.OptNoTrans
//   v := alpha * A' * u + beta * v  if trans is linalg.OptTrans
//
type MatrixA interface {
	Af(u, v *matrix.FloatMatrix, alpha, beta float64, trans la.Option) error
}

//
// MatrixVarA provides interface for extended customization with non-matrix type
// primal and dual variables.
//
type MatrixVarA interface {
	Af(u, v MatrixVariable, alpha, beta float64, trans la.Option) error
}

// Public interface to provide custom P matrix-vector products
//
// The call Pf(u, v, alpha, beta) should evaluate the matrix-vectors product.
//
//   v := alpha * P * u + beta * v.
//
type MatrixP interface {
	Pf(u, v *matrix.FloatMatrix, alpha, beta float64) error
}

//
// MatrixVarP provides interface for extended customization with non-matrix type
// primal and dual variables.
//
type MatrixVarP interface {
	Pf(u, v MatrixVariable, alpha, beta float64) error
}

// ConvexProg is an interface that handles the following functions.
//
// F0() returns a tuple (mnl, x0, error).
//
//  - mnl is the number of nonlinear inequality constraints.
//  - x0 is a point in the domain of f.
//
// F1(x) returns a tuple (f, Df, error).
//
//  - f is a matrix of size (mnl, 1) containing f(x).
//  - Df is a matrix of size (mnl, n), containing the derivatives of f at x.
//    Df[k,:] is the transpose of the gradient of fk at x.
//    If x is not in dom f, F1(x) returns (nil, nil, error)
//
// F2(x, z) with z a positive matrix of size (mnl,1), returns a tuple (f, Df, H, error).
//
//   - f and Df are defined as above.
//   - H is a matrix of size (n,n).  The lower triangular part of H contains the
//     lower triangular part of sum_k z[k] * Hk where Hk is the Hessian of fk at x.
//
// When F2() is called, it can be assumed that x is dom f.
//
//
type ConvexProg interface {
	// Returns (mnl, x0) where mln number of nonlinear inequality constraints
	// and x0 is a point in the domain of f.
	F0() (mnl int, x0 *matrix.FloatMatrix, err error)

	// Returns a tuple (f, Df) where f is of size (mnl, 1) containing f(x)
	// Df is matrix of size (mnl, n) containing the derivatives of f at x:
	// Df[k,:] is the transpose of the gradient of fk at x. If x is not in
	// domf, return non-nil error.
	F1(x *matrix.FloatMatrix) (f, Df *matrix.FloatMatrix, err error)

	// F(x, z) with z a positive  matrix of size (mnl, 1). Return a tuple
	// (f, Df, H), where f, Df as above. H is matrix of size (n, n).
	F2(x, z *matrix.FloatMatrix) (f, Df, H *matrix.FloatMatrix, err error)
}

// ConvexVarProg provides interface for extended customization with non-matrix type
// primal and dual variables.
//
type ConvexVarProg interface {
	// Returns (mnl, x0) where mln number of nonlinear inequality constraints
	// and x0 is a point in the domain of f.
	F0() (mnl int, x0 MatrixVariable, err error)

	// Returns a tuple (f, Df) where f contains f(x) and Df is instance of
	// MatrixVarDf interface.
	F1(x MatrixVariable) (f MatrixVariable, Df MatrixVarDf, err error)

	// F(x, z) with z a positive matrix variable. Return a tuple
	// (f, Df, H), where f, Df as above. H is an instance of MatrixVarH interface.
	F2(x MatrixVariable, z *matrix.FloatMatrix) (f MatrixVariable, Df MatrixVarDf, H MatrixVarH, err error)
}

// Public interface to provide custom Df matrix-vector products
//
// The call Df(u, v, alpha, beta, trans) should evaluate the matrix-vectors products
//
//   v := alpha * Df(x) * u + beta * v if trans is linalg.OptNoTrans
//   v := alpha * Df(x)' * u + beta * v if trans is linalg.OptTrans.
//
type MatrixDf interface {
	Df(u, v *matrix.FloatMatrix, alpha, beta float64, trans la.Option) error
}

// MatrixVarDf provides interface for extended customization with non-matrix type
// primal and dual variables.
//
type MatrixVarDf interface {
	Df(u, v MatrixVariable, alpha, beta float64, trans la.Option) error
}

// Public interface to provide custom H matrix-vector products
//
// The call H(u, v, alpha, beta) should evaluate the matrix-vectors product
//
//    v := alpha * H * u + beta * v.
//
type MatrixH interface {
	Hf(u, v *matrix.FloatMatrix, alpha, beta float64) error
}

// MatrixVarH provides interface for extended customization with non-matrix type
// primal and dual variables.
//
type MatrixVarH interface {
	Hf(u, v MatrixVariable, alpha, beta float64) error
}

// MatrixVariable interface for any type used to represent primal variables
// and the dual variables as something else than one-column float matrices.
//
type MatrixVariable interface {
	// Provide internal matrix value
	Matrix() *matrix.FloatMatrix
	// Create a new copy
	Copy() MatrixVariable
	// Computes v := alpha*u + v for a scalar alpha and vectors u and v.
	Axpy(v MatrixVariable, alpha float64)
	// Return the inner product of two vectors u and v in a vector space.
	Dot(v MatrixVariable) float64
	// Computes u := alpha*u for a scalar alpha and vectors u in a vector space.
	Scal(alpha float64)
	// Implement checkpnt.Verifiable to allow checkpoint checking
	Verify(dataline string) float64
	// Implement checkpnt.Verifiable to allow checkpoint checking
	ShowError(dataline string)
	// Convert to string for printing
	String() string
}

// Copies x to y.
func mCopy(x, y MatrixVariable) {
	y.Scal(0.0)
	x.Axpy(y, 1.0)
}

// Implementations of interfaces for standard matrix valued parameters

// Impelements MatrixVariable interface standard matrix valued variable
type matrixVar struct {
	val *matrix.FloatMatrix
}

// Internal matrix
func (u *matrixVar) Matrix() *matrix.FloatMatrix {
	return u.val
}

// New copy
func (u *matrixVar) Copy() MatrixVariable {
	return &matrixVar{u.val.Copy()}
}

// Inner product
func (u *matrixVar) Dot(v MatrixVariable) float64 {
	if y, ok := v.(*matrixVar); ok {
		return blas.DotFloat(u.val, y.val)
	}
	return 0.0
}

func (u *matrixVar) Axpy(v MatrixVariable, alpha float64) {
	if y, ok := v.(*matrixVar); ok {
		blas.AxpyFloat(u.val, y.val, alpha)
	}
	return
}

func (u *matrixVar) Scal(alpha float64) {
	blas.ScalFloat(u.val, alpha)
}

func (u *matrixVar) String() string {
	return u.val.ToString("%.7f")
}

func (u *matrixVar) Verify(dataline string) float64 {
	diff := 0.0
	refval, _ := matrix.FloatParse(dataline)
	diff += blas.Nrm2Float(matrix.Minus(u.val, refval))
	return diff
}

func (u *matrixVar) ShowError(dataline string) {
	refval, _ := matrix.FloatParse(dataline)
	em := matrix.Minus(u.val, refval)
	r, _ := matrix.FloatMatrixStacked(matrix.StackRight, u.val, refval, em)
	fmt.Printf("my data | ref.data | diff \n%v\n", r.ToString("%.4e"))
}

// Implements MatrixVarG interface for matrix valued G
type matrixVarG struct {
	mG   *matrix.FloatMatrix
	dims *sets.DimensionSet
}

func (g *matrixVarG) Gf(u, v MatrixVariable, alpha, beta float64, trans la.Option) (err error) {
	err = sgemv(g.mG, u.Matrix(), v.Matrix(), alpha, beta, g.dims, trans)
	return
}

// Implements MatrixVarG interface for MatrixG interface and matrix valued u, v
type matrixIfG struct {
	mG MatrixG
}

func (g *matrixIfG) Gf(u, v MatrixVariable, alpha, beta float64, trans la.Option) (err error) {
	err = g.mG.Gf(u.Matrix(), v.Matrix(), alpha, beta, trans)
	return
}

// Implements MatrixVarA interface for matrix valued A, u, v
type matrixVarA struct {
	mA *matrix.FloatMatrix
}

func (a *matrixVarA) Af(u, v MatrixVariable, alpha, beta float64, trans la.Option) (err error) {
	err = blas.GemvFloat(a.mA, u.Matrix(), v.Matrix(), alpha, beta, trans)
	return
}

// Implements MatrixVarA interface for MatrixA interface
type matrixIfA struct {
	mA MatrixA
}

func (a *matrixIfA) Af(u, v MatrixVariable, alpha, beta float64, trans la.Option) (err error) {
	err = a.mA.Af(u.Matrix(), v.Matrix(), alpha, beta, trans)
	return
}

// Implements MatrixVarA interface for matrix valued A, u, v
type matrixVarP struct {
	mP *matrix.FloatMatrix
}

func (p *matrixVarP) Pf(u, v MatrixVariable, alpha, beta float64) (err error) {
	err = blas.SymvFloat(p.mP, u.Matrix(), v.Matrix(), alpha, beta)
	return
}

// Implements MatrixVarA interface for MatrixA interface
type matrixIfP struct {
	mP MatrixP
}

func (p *matrixIfP) Pf(u, v MatrixVariable, alpha, beta float64) (err error) {
	err = p.mP.Pf(u.Matrix(), v.Matrix(), alpha, beta)
	return
}

// Implements ConvexVarProg interface for standard matrices
type convexVarProg struct {
	cp ConvexProg
}

func (p *convexVarProg) F0() (mnl int, x0 MatrixVariable, err error) {
	x0 = nil
	var mx0 *matrix.FloatMatrix
	mnl, mx0, err = p.cp.F0()
	if err != nil {
		return
	}
	x0 = &matrixVar{mx0}
	return
}

func (p *convexVarProg) F1(x MatrixVariable) (f MatrixVariable, Df MatrixVarDf, err error) {
	f = nil
	Df = nil
	var f0, Df0 *matrix.FloatMatrix
	f0, Df0, err = p.cp.F1(x.Matrix())
	if err != nil {
		return
	}
	f = &matrixVar{f0}
	Df = &matrixVarDf{Df0}
	return
}

func (p *convexVarProg) F2(x MatrixVariable, z *matrix.FloatMatrix) (f MatrixVariable, Df MatrixVarDf, H MatrixVarH, err error) {
	f = nil
	Df = nil
	H = nil
	var f0, Df0, H0 *matrix.FloatMatrix
	f0, Df0, H0, err = p.cp.F2(x.Matrix(), z)
	if err != nil {
		return
	}
	f = &matrixVar{f0}
	Df = &matrixVarDf{Df0}
	H = &matrixVarH{H0}
	return
}

// Implements MatrixVarDf interface for standard matrix
type matrixVarDf struct {
	df *matrix.FloatMatrix
}

func (d *matrixVarDf) Df(u, v MatrixVariable, alpha, beta float64, trans la.Option) error {
	return blas.GemvFloat(d.df, u.Matrix(), v.Matrix(), alpha, beta, trans)
}

// Implements MatrixVarH interface for standard matrix
type matrixVarH struct {
	h *matrix.FloatMatrix
}

func (hf *matrixVarH) Hf(u, v MatrixVariable, alpha, beta float64) error {
	return blas.SymvFloat(hf.h, u.Matrix(), v.Matrix(), alpha, beta)
}

// Local Variables:
// tab-width: 4
// End:
