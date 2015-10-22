// Copyright (c) Harri Rautila, 2012

// This file is part of github.com/hrautila/cvx package.
// It is free software, distributed under the terms of GNU Lesser General Public
// License Version 3, or any later version. See the COPYING tile included in this archive.

package cvx

import (
	"errors"
	"fmt"
	"github.com/henrylee2cn/algorithm/cvx/sets"
	"github.com/henrylee2cn/algorithm/matrix"
	"math"
)

// Solves a pair of primal and dual LPs
//
//      minimize    c'*x
//      subject to  G*x + s = h
//                  A*x = b
//                  s >= 0
//
//      maximize    -h'*z - b'*y
//      subject to  G'*z + A'*y + c = 0
//                    z >= 0.
//
func Lp(c, G, h, A, b *matrix.FloatMatrix, solopts *SolverOptions,
	primalstart, dualstart *sets.FloatMatrixSet) (sol *Solution, err error) {

	if c == nil {
		err = errors.New("'c' must a column matrix")
		return
	}
	n := c.Rows()
	if n < 1 {
		err = errors.New("Number of variables must be at least 1")
		return
	}
	if G == nil || G.Cols() != n {
		err = errors.New(fmt.Sprintf("'G' must be matrix with %d columns", n))
		return
	}
	m := G.Rows()
	if h == nil || !h.SizeMatch(m, 1) {
		err = errors.New(fmt.Sprintf("'h' must be matrix of size (%d,1)", m))
		return
	}
	if A == nil {
		A = matrix.FloatZeros(0, n)
	}
	if A.Cols() != n {
		err = errors.New(fmt.Sprintf("'A' must be matrix with %d columns", n))
		return
	}
	p := A.Rows()
	if b == nil {
		b = matrix.FloatZeros(0, 1)
	}
	if !b.SizeMatch(p, 1) {
		err = errors.New(fmt.Sprintf("'b' must be matrix of size (%d,1)", p))
		return
	}
	dims := sets.NewDimensionSet("l", "q", "s")
	dims.Set("l", []int{m})

	return ConeLp(c, G, h, A, b, dims, solopts, primalstart, dualstart)
}

// Solves a quadratic program
//
//        minimize    (1/2)*x'*P*x + q'*x
//        subject to  G*x <= h
//                    A*x = b.
//
//
func Qp(P, q, G, h, A, b *matrix.FloatMatrix, solopts *SolverOptions,
	initvals *sets.FloatMatrixSet) (sol *Solution, err error) {

	sol = nil
	if P == nil || P.Rows() != P.Cols() {
		err = errors.New("'P' must a non-nil square matrix")
		return
	}
	if q == nil {
		err = errors.New("'q' must a non-nil matrix")
		return
	}
	if q.Rows() != P.Rows() || q.Cols() > 1 {
		err = errors.New(fmt.Sprintf("'q' must be matrix of size (%d,1)", P.Rows()))
		return
	}
	if G == nil {
		G = matrix.FloatZeros(0, P.Rows())
	}
	if G.Cols() != P.Rows() {
		err = errors.New(fmt.Sprintf("'G' must be matrix of %d columns", P.Rows()))
		return
	}
	if h == nil {
		h = matrix.FloatZeros(G.Rows(), 1)
	}
	if h.Rows() != G.Rows() || h.Cols() > 1 {
		err = errors.New(fmt.Sprintf("'h' must be matrix of size (%d,1)", G.Rows()))
		return
	}
	if A == nil {
		A = matrix.FloatZeros(0, P.Rows())
	}
	if A.Cols() != P.Rows() {
		err = errors.New(fmt.Sprintf("'A' must be matrix of %d columns", P.Rows()))
		return
	}
	if b == nil {
		b = matrix.FloatZeros(A.Rows(), 1)
	}
	if b.Rows() != A.Rows() {
		err = errors.New(fmt.Sprintf("'b' must be matrix of size (%d,1)", A.Rows()))
		return
	}
	return ConeQp(P, q, G, h, A, b, nil, solopts, initvals)
}

// Solves a pair of primal and dual SOCPs
//
//     minimize    c'*x
//     subject to  Gl*x + sl = hl
//                 Gq[k]*x + sq[k] = hq[k],  k = 0, ..., N-1
//                 A*x = b
//                 sl >= 0,
//                 sq[k] >= 0, k = 0, ..., N-1
//
//     maximize   -hl'*z - sum_k hq[k]'*zq[k] - b'*y
//     subject to  Gl'*zl + sum_k Gq[k]'*zq[k] + A'*y + c = 0
//                 zl >= 0,  zq[k] >= 0, k = 0, ..., N-1.
//
// The inequalities sl >= 0 and zl >= 0 are elementwise vector
// inequalities.  The inequalities sq[k] >= 0, zq[k] >= 0 are second
// order cone inequalities, i.e., equivalent to
//
//     sq[k][0] >= || sq[k][1:] ||_2,  zq[k][0] >= || zq[k][1:] ||_2.
//
func Socp(c, Gl, hl, A, b *matrix.FloatMatrix, Ghq *sets.FloatMatrixSet, solopts *SolverOptions,
	primalstart, dualstart *sets.FloatMatrixSet) (sol *Solution, err error) {
	if c == nil {
		err = errors.New("'c' must a column matrix")
		return
	}
	n := c.Rows()
	if n < 1 {
		err = errors.New("Number of variables must be at least 1")
		return
	}
	if Gl == nil {
		Gl = matrix.FloatZeros(0, n)
	}
	if Gl.Cols() != n {
		err = errors.New(fmt.Sprintf("'G' must be matrix with %d columns", n))
		return
	}
	ml := Gl.Rows()
	if hl == nil {
		hl = matrix.FloatZeros(0, 1)
	}
	if !hl.SizeMatch(ml, 1) {
		err = errors.New(fmt.Sprintf("'hl' must be matrix of size (%d,1)", ml))
		return
	}
	Gqset := Ghq.At("Gq")
	mq := make([]int, 0)
	for i, Gq := range Gqset {
		if Gq.Cols() != n {
			err = errors.New(fmt.Sprintf("'Gq' must be list of matrices with %d columns", n))
			return
		}
		if Gq.Rows() == 0 {
			err = errors.New(fmt.Sprintf("the number of rows of 'Gq[%d]' is zero", i))
			return
		}
		mq = append(mq, Gq.Rows())
	}
	hqset := Ghq.At("hq")
	if len(Gqset) != len(hqset) {
		err = errors.New(fmt.Sprintf("'hq' must be a list of %d matrices", len(Gqset)))
		return
	}
	for i, hq := range hqset {
		if !hq.SizeMatch(Gqset[i].Rows(), 1) {
			s := fmt.Sprintf("hq[%d] has size (%d,%d). Expected size is (%d,1)",
				i, hq.Rows(), hq.Cols(), Gqset[i].Rows())
			err = errors.New(s)
			return
		}
	}
	if A == nil {
		A = matrix.FloatZeros(0, n)
	}
	if A.Cols() != n {
		err = errors.New(fmt.Sprintf("'A' must be matrix with %d columns", n))
		return
	}
	p := A.Rows()
	if b == nil {
		b = matrix.FloatZeros(0, 1)
	}
	if !b.SizeMatch(p, 1) {
		err = errors.New(fmt.Sprintf("'b' must be matrix of size (%d,1)", p))
		return
	}
	dims := sets.NewDimensionSet("l", "q", "s")
	dims.Set("l", []int{ml})
	dims.Set("q", mq)
	//N := dims.Sum("l", "q")

	hargs := make([]*matrix.FloatMatrix, 0, len(hqset)+1)
	hargs = append(hargs, hl)
	hargs = append(hargs, hqset...)
	h, indh := matrix.FloatMatrixStacked(matrix.StackDown, hargs...)

	Gargs := make([]*matrix.FloatMatrix, 0, len(Gqset)+1)
	Gargs = append(Gargs, Gl)
	Gargs = append(Gargs, Gqset...)
	G, indg := matrix.FloatMatrixStacked(matrix.StackDown, Gargs...)

	var pstart, dstart *sets.FloatMatrixSet = nil, nil
	if primalstart != nil {
		pstart = sets.NewFloatSet("x", "s")
		pstart.Set("x", primalstart.At("x")[0])
		slset := primalstart.At("sl")
		margs := make([]*matrix.FloatMatrix, 0, len(slset)+1)
		margs = append(margs, primalstart.At("s")[0])
		margs = append(margs, slset...)
		sl, _ := matrix.FloatMatrixStacked(matrix.StackDown, margs...)
		pstart.Set("s", sl)
	}

	if dualstart != nil {
		dstart = sets.NewFloatSet("y", "z")
		dstart.Set("y", dualstart.At("y")[0])
		zlset := primalstart.At("zl")
		margs := make([]*matrix.FloatMatrix, 0, len(zlset)+1)
		margs = append(margs, dualstart.At("z")[0])
		margs = append(margs, zlset...)
		zl, _ := matrix.FloatMatrixStacked(matrix.StackDown, margs...)
		dstart.Set("z", zl)
	}

	sol, err = ConeLp(c, G, h, A, b, dims, solopts, pstart, dstart)
	// unpack sol.Result
	if err == nil {
		s := sol.Result.At("s")[0]
		sl := matrix.FloatVector(s.FloatArray()[:ml])
		sol.Result.Append("sl", sl)
		ind := ml
		for _, k := range indh[1:] {
			sk := matrix.FloatVector(s.FloatArray()[ind : ind+k])
			sol.Result.Append("sq", sk)
			ind += k
		}

		z := sol.Result.At("z")[0]
		zl := matrix.FloatVector(z.FloatArray()[:ml])
		sol.Result.Append("zl", zl)
		ind = ml
		for _, k := range indg[1:] {
			zk := matrix.FloatVector(z.FloatArray()[ind : ind+k])
			sol.Result.Append("zq", zk)
			ind += k
		}
	}
	sol.Result.Remove("s")
	sol.Result.Remove("z")

	return
}

// Solves a pair of primal and dual SDPs
//
//        minimize    c'*x
//        subject to  Gl*x + sl = hl
//                    mat(Gs[k]*x) + ss[k] = hs[k], k = 0, ..., N-1
//                    A*x = b
//                    sl >= 0,  ss[k] >= 0, k = 0, ..., N-1
//
//        maximize    -hl'*z - sum_k trace(hs[k]*zs[k]) - b'*y
//        subject to  Gl'*zl + sum_k Gs[k]'*vec(zs[k]) + A'*y + c = 0
//                    zl >= 0,  zs[k] >= 0, k = 0, ..., N-1.
//
// The inequalities sl >= 0 and zl >= 0 are elementwise vector
// inequalities.  The inequalities ss[k] >= 0, zs[k] >= 0 are matrix
// inequalities, i.e., the symmetric matrices ss[k] and zs[k] must be
// positive semidefinite.  mat(Gs[k]*x) is the symmetric matrix X with
// X[:] = Gs[k]*x.  For a symmetric matrix, zs[k], vec(zs[k]) is the
// vector zs[k][:].
//
func Sdp(c, Gl, hl, A, b *matrix.FloatMatrix, Ghs *sets.FloatMatrixSet, solopts *SolverOptions,
	primalstart, dualstart *sets.FloatMatrixSet) (sol *Solution, err error) {
	if c == nil {
		err = errors.New("'c' must a column matrix")
		return
	}
	n := c.Rows()
	if n < 1 {
		err = errors.New("Number of variables must be at least 1")
		return
	}
	if Gl == nil {
		Gl = matrix.FloatZeros(0, n)
	}
	if Gl.Cols() != n {
		err = errors.New(fmt.Sprintf("'G' must be matrix with %d columns", n))
		return
	}
	ml := Gl.Rows()
	if hl == nil {
		hl = matrix.FloatZeros(0, 1)
	}
	if !hl.SizeMatch(ml, 1) {
		err = errors.New(fmt.Sprintf("'hl' must be matrix of size (%d,1)", ml))
		return
	}
	Gsset := Ghs.At("Gs")
	ms := make([]int, 0)
	for i, Gs := range Gsset {
		if Gs.Cols() != n {
			err = errors.New(fmt.Sprintf("'Gs' must be list of matrices with %d columns", n))
			return
		}
		sz := int(math.Sqrt(float64(Gs.Rows())))
		if Gs.Rows() != sz*sz {
			err = errors.New(fmt.Sprintf("the squareroot of the number of rows of 'Gq[%d]' is not an integer", i))
			return
		}
		ms = append(ms, sz)
	}

	hsset := Ghs.At("hs")
	if len(Gsset) != len(hsset) {
		err = errors.New(fmt.Sprintf("'hs' must be a list of %d matrices", len(Gsset)))
		return
	}
	for i, hs := range hsset {
		if !hs.SizeMatch(ms[i], ms[i]) {
			s := fmt.Sprintf("hq[%d] has size (%d,%d). Expected size is (%d,%d)",
				i, hs.Rows(), hs.Cols(), ms[i], ms[i])
			err = errors.New(s)
			return
		}
	}
	if A == nil {
		A = matrix.FloatZeros(0, n)
	}
	if A.Cols() != n {
		err = errors.New(fmt.Sprintf("'A' must be matrix with %d columns", n))
		return
	}
	p := A.Rows()
	if b == nil {
		b = matrix.FloatZeros(0, 1)
	}
	if !b.SizeMatch(p, 1) {
		err = errors.New(fmt.Sprintf("'b' must be matrix of size (%d,1)", p))
		return
	}
	dims := sets.NewDimensionSet("l", "q", "s")
	dims.Set("l", []int{ml})
	dims.Set("s", ms)
	N := dims.Sum("l") + dims.SumSquared("s")

	// Map hs matrices to h vector
	h := matrix.FloatZeros(N, 1)
	h.SetIndexesFromArray(hl.FloatArray()[:ml], matrix.MakeIndexSet(0, ml, 1)...)
	ind := ml
	for k, hs := range hsset {
		h.SetIndexesFromArray(hs.FloatArray(), matrix.MakeIndexSet(ind, ind+ms[k]*ms[k], 1)...)
		ind += ms[k] * ms[k]
	}

	Gargs := make([]*matrix.FloatMatrix, 0)
	Gargs = append(Gargs, Gl)
	Gargs = append(Gargs, Gsset...)
	G, sizeg := matrix.FloatMatrixStacked(matrix.StackDown, Gargs...)

	var pstart, dstart *sets.FloatMatrixSet = nil, nil
	if primalstart != nil {
		pstart = sets.NewFloatSet("x", "s")
		pstart.Set("x", primalstart.At("x")[0])
		slset := primalstart.At("sl")
		margs := make([]*matrix.FloatMatrix, 0, len(slset)+1)
		margs = append(margs, primalstart.At("s")[0])
		margs = append(margs, slset...)
		sl, _ := matrix.FloatMatrixStacked(matrix.StackDown, margs...)
		pstart.Set("s", sl)
	}

	if dualstart != nil {
		dstart = sets.NewFloatSet("y", "z")
		dstart.Set("y", dualstart.At("y")[0])
		zlset := primalstart.At("zl")
		margs := make([]*matrix.FloatMatrix, 0, len(zlset)+1)
		margs = append(margs, dualstart.At("z")[0])
		margs = append(margs, zlset...)
		zl, _ := matrix.FloatMatrixStacked(matrix.StackDown, margs...)
		dstart.Set("z", zl)
	}

	//fmt.Printf("h=\n%v\n", h.ToString("%.3f"))
	//fmt.Printf("G=\n%v\n", G.ToString("%.3f"))

	sol, err = ConeLp(c, G, h, A, b, dims, solopts, pstart, dstart)
	// unpack sol.Result
	if err == nil {
		s := sol.Result.At("s")[0]
		sl := matrix.FloatVector(s.FloatArray()[:ml])
		sol.Result.Append("sl", sl)
		ind := ml
		for _, m := range ms {
			sk := matrix.FloatNew(m, m, s.FloatArray()[ind:ind+m*m])
			sol.Result.Append("ss", sk)
			ind += m * m
		}

		z := sol.Result.At("z")[0]
		zl := matrix.FloatVector(s.FloatArray()[:ml])
		sol.Result.Append("zl", zl)
		ind = ml
		for i, k := range sizeg[1:] {
			zk := matrix.FloatNew(ms[i], ms[i], z.FloatArray()[ind:ind+k])
			sol.Result.Append("zs", zk)
			ind += k
		}
	}
	sol.Result.Remove("s")
	sol.Result.Remove("z")

	return

}

// Local Variables:
// tab-width: 4
// End:
