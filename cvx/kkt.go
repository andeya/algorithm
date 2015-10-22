// Copyright (c) Harri Rautila, 2012

// This file is part of github.com/hrautila/cvx package.
// It is free software, distributed under the terms of GNU Lesser General Public
// License Version 3, or any later version. See the COPYING tile included in this archive.

package cvx

import (
	"errors"
	"github.com/henrylee2cn/algorithm/cvx/checkpnt"
	"github.com/henrylee2cn/algorithm/cvx/sets"
	la "github.com/henrylee2cn/algorithm/linalg"
	"github.com/henrylee2cn/algorithm/linalg/blas"
	"github.com/henrylee2cn/algorithm/linalg/lapack"
	"github.com/henrylee2cn/algorithm/matrix"
	//"fmt"
	//"math"
)

func setDiagonal(M *matrix.FloatMatrix, srow, scol, erow, ecol int, val float64) {
	for i := srow; i < erow; i++ {
		if i < ecol {
			M.SetAt(i, i, val)
		}
	}
}

// Solution of KKT equations by a dense LDL factorization of the
// 3 x 3 system.
//
// Returns a function that (1) computes the LDL factorization of
//
// [ H           A'   GG'*W^{-1} ]
// [ A           0    0          ],
// [ W^{-T}*GG   0   -I          ]
//
// given H, Df, W, where GG = [Df; G], and (2) returns a function for
// solving
//
//  [ H     A'   GG'   ]   [ ux ]   [ bx ]
//  [ A     0    0     ] * [ uy ] = [ by ].
//  [ GG    0   -W'*W  ]   [ uz ]   [ bz ]
//
// H is n x n,  A is p x n, Df is mnl x n, G is N x n where
// N = dims['l'] + sum(dims['q']) + sum( k**2 for k in dims['s'] ).
//
func kktLdl(G *matrix.FloatMatrix, dims *sets.DimensionSet, A *matrix.FloatMatrix, mnl int) (kktFactor, error) {

	p, n := A.Size()
	ldK := n + p + mnl + dims.At("l")[0] + dims.Sum("q") + dims.SumPacked("s")
	K := matrix.FloatZeros(ldK, ldK)
	ipiv := make([]int32, ldK)
	u := matrix.FloatZeros(ldK, 1)
	g := matrix.FloatZeros(mnl+G.Rows(), 1)
	//checkpnt.AddMatrixVar("u", u)
	//checkpnt.AddMatrixVar("K", K)

	factor := func(W *sets.FloatMatrixSet, H, Df *matrix.FloatMatrix) (KKTFunc, error) {
		var err error = nil
		// Zero K for each call.
		blas.ScalFloat(K, 0.0)
		if H != nil {
			K.SetSubMatrix(0, 0, H)
		}
		K.SetSubMatrix(n, 0, A)
		for k := 0; k < n; k++ {
			// g is (mnl + G.Rows(), 1) matrix, Df is (mnl, n), G is (N, n)
			if mnl > 0 {
				// set values g[0:mnl] = Df[,k]
				g.SetIndexesFromArray(Df.GetColumnArray(k, nil), matrix.MakeIndexSet(0, mnl, 1)...)
			}
			// set values g[mnl:] = G[,k]
			g.SetIndexesFromArray(G.GetColumnArray(k, nil), matrix.MakeIndexSet(mnl, mnl+g.Rows(), 1)...)
			scale(g, W, true, true)
			if err != nil {
				//fmt.Printf("scale error: %s\n", err)
			}
			pack(g, K, dims, &la.IOpt{"mnl", mnl}, &la.IOpt{"offsety", k*ldK + n + p})
		}
		setDiagonal(K, n+p, n+n, ldK, ldK, -1.0)
		err = lapack.Sytrf(K, ipiv)
		if err != nil {
			return nil, err
		}

		solve := func(x, y, z *matrix.FloatMatrix) (err error) {
			// Solve
			//
			//     [ H          A'   GG'*W^{-1} ]   [ ux   ]   [ bx        ]
			//     [ A          0    0          ] * [ uy   [ = [ by        ]
			//     [ W^{-T}*GG  0   -I          ]   [ W*uz ]   [ W^{-T}*bz ]
			//
			// and return ux, uy, W*uz.
			//
			// On entry, x, y, z contain bx, by, bz.  On exit, they contain
			// the solution ux, uy, W*uz.
			err = nil
			blas.Copy(x, u)
			blas.Copy(y, u, &la.IOpt{"offsety", n})
			err = scale(z, W, true, true)
			if err != nil {
				return
			}
			err = pack(z, u, dims, &la.IOpt{"mnl", mnl}, &la.IOpt{"offsety", n + p})
			if err != nil {
				return
			}

			err = lapack.Sytrs(K, u, ipiv)
			if err != nil {
				return
			}

			blas.Copy(u, x, &la.IOpt{"n", n})
			blas.Copy(u, y, &la.IOpt{"n", p}, &la.IOpt{"offsetx", n})
			err = unpack(u, z, dims, &la.IOpt{"mnl", mnl}, &la.IOpt{"offsetx", n + p})
			return
		}
		return solve, err
	}
	return factor, nil
}

// Solution of KKT equations with zero 1,1 block, by eliminating the
// equality constraints via a QR factorization, and solving the
// reduced KKT system by another QR factorization.
//
// Computes the QR factorization
//
//        A' = [Q1, Q2] * [R1; 0]
//
// and returns a function that (1) computes the QR factorization
//
//        W^{-T} * G * Q2 = Q3 * R3
//
// (with columns of W^{-T}*G in packed storage), and (2) returns a function for solving
//
//        [ 0    A'   G'    ]   [ ux ]   [ bx ]
//        [ A    0    0     ] * [ uy ] = [ by ].
//        [ G    0   -W'*W  ]   [ uz ]   [ bz ]
//
// A is p x n and G is N x n where N = dims['l'] + sum(dims['q']) +
// sum( k**2 for k in dims['s'] ).
//
func kktQr(G *matrix.FloatMatrix, dims *sets.DimensionSet, A *matrix.FloatMatrix, mnl int) (kktFactor, error) {

	p, n := A.Size()
	cdim := dims.Sum("l", "q") + dims.SumSquared("s")
	cdim_pckd := dims.Sum("l", "q") + dims.SumPacked("s")

	QA := A.Transpose()
	tauA := matrix.FloatZeros(p, 1)
	lapack.Geqrf(QA, tauA)

	Gs := matrix.FloatZeros(cdim, n)
	tauG := matrix.FloatZeros(n-p, 1)
	u := matrix.FloatZeros(cdim_pckd, 1)
	vv := matrix.FloatZeros(n, 1)
	w := matrix.FloatZeros(cdim_pckd, 1)
	checkpnt.AddMatrixVar("tauA", tauA)
	checkpnt.AddMatrixVar("tauG", tauG)
	checkpnt.AddMatrixVar("Gs", Gs)
	checkpnt.AddMatrixVar("qr_u", u)
	checkpnt.AddMatrixVar("qr_vv", vv)

	factor := func(W *sets.FloatMatrixSet, H, Df *matrix.FloatMatrix) (KKTFunc, error) {
		var err error = nil
		minor := 0
		if !checkpnt.MinorEmpty() {
			minor = checkpnt.MinorTop()
		}

		// Gs = W^{-T}*G, in packed storage.
		blas.Copy(G, Gs)
		//checkpnt.Check("00factor_qr", minor)
		scale(Gs, W, true, true)
		//checkpnt.Check("01factor_qr", minor)
		pack2(Gs, dims, 0)
		//checkpnt.Check("02factor_qr", minor)

		// Gs := [ Gs1, Gs2 ]
		//     = Gs * [ Q1, Q2 ]
		lapack.Ormqr(QA, tauA, Gs, la.OptRight, &la.IOpt{"m", cdim_pckd})
		//checkpnt.Check("03factor_qr", minor)

		// QR factorization Gs2 := [ Q3, Q4 ] * [ R3; 0 ]
		lapack.Geqrf(Gs, tauG, &la.IOpt{"n", n - p}, &la.IOpt{"m", cdim_pckd},
			&la.IOpt{"offseta", Gs.Rows() * p})
		checkpnt.Check("10factor_qr", minor)

		solve := func(x, y, z *matrix.FloatMatrix) (err error) {
			// On entry, x, y, z contain bx, by, bz.  On exit, they
			// contain the solution x, y, W*z of
			//
			//     [ 0         A'  G'*W^{-1} ]   [ x   ]   [bx       ]
			//     [ A         0   0         ] * [ y   ] = [by       ].
			//     [ W^{-T}*G  0   -I        ]   [ W*z ]   [W^{-T}*bz]
			//
			// The system is solved in five steps:
			//
			//       w := W^{-T}*bz - Gs1*R1^{-T}*by
			//       u := R3^{-T}*Q2'*bx + Q3'*w
			//     W*z := Q3*u - w
			//       y := R1^{-1} * (Q1'*bx - Gs1'*(W*z))
			//       x := [ Q1, Q2 ] * [ R1^{-T}*by;  R3^{-1}*u ]

			minor := 0
			if !checkpnt.MinorEmpty() {
				minor = checkpnt.MinorTop()
			}

			// w := W^{-T} * bz in packed storage
			scale(z, W, true, true)
			pack(z, w, dims)
			//checkpnt.Check("00solve_qr", minor)

			// vv := [ Q1'*bx;  R3^{-T}*Q2'*bx ]
			blas.Copy(x, vv)
			lapack.Ormqr(QA, tauA, vv, la.OptTrans)
			lapack.Trtrs(Gs, vv, la.OptUpper, la.OptTrans, &la.IOpt{"n", n - p},
				&la.IOpt{"offseta", Gs.Rows() * p}, &la.IOpt{"offsetb", p})
			//checkpnt.Check("10solve_qr", minor)

			// x[:p] := R1^{-T} * by
			blas.Copy(y, x)
			lapack.Trtrs(QA, x, la.OptUpper, la.OptTrans, &la.IOpt{"n", p})
			//checkpnt.Check("20solve_qr", minor)

			// w := w - Gs1 * x[:p]
			//    = W^{-T}*bz - Gs1*by
			blas.GemvFloat(Gs, x, w, -1.0, 1.0, &la.IOpt{"n", p}, &la.IOpt{"m", cdim_pckd})
			//checkpnt.Check("30solve_qr", minor)

			// u := [ Q3'*w + v[p:];  0 ]
			//    = [ Q3'*w + R3^{-T}*Q2'*bx; 0 ]
			blas.Copy(w, u)
			lapack.Ormqr(Gs, tauG, u, la.OptTrans, &la.IOpt{"k", n - p},
				&la.IOpt{"offseta", Gs.Rows() * p}, &la.IOpt{"m", cdim_pckd})
			blas.AxpyFloat(vv, u, 1.0, &la.IOpt{"offsetx", p}, &la.IOpt{"n", n - p})
			blas.ScalFloat(u, 0.0, &la.IOpt{"offset", n - p})
			//checkpnt.Check("40solve_qr", minor)

			// x[p:] := R3^{-1} * u[:n-p]
			blas.Copy(u, x, &la.IOpt{"offsety", p}, &la.IOpt{"n", n - p})
			lapack.Trtrs(Gs, x, la.OptUpper, &la.IOpt{"n", n - p},
				&la.IOpt{"offset", Gs.Rows() * p}, &la.IOpt{"offsetb", p})
			//checkpnt.Check("50solve_qr", minor)

			// x is now [ R1^{-T}*by;  R3^{-1}*u[:n-p] ]
			// x := [Q1 Q2]*x
			lapack.Ormqr(QA, tauA, x)
			//checkpnt.Check("60solve_qr", minor)

			// u := [Q3, Q4] * u - w
			//    = Q3 * u[:n-p] - w
			lapack.Ormqr(Gs, tauG, u, &la.IOpt{"k", n - p}, &la.IOpt{"m", cdim_pckd},
				&la.IOpt{"offseta", Gs.Rows() * p})
			blas.AxpyFloat(w, u, -1.0)
			//checkpnt.Check("70solve_qr", minor)

			// y := R1^{-1} * ( v[:p] - Gs1'*u )
			//    = R1^{-1} * ( Q1'*bx - Gs1'*u )
			blas.Copy(vv, y, &la.IOpt{"n", p})
			blas.GemvFloat(Gs, u, y, -1.0, 1.0, &la.IOpt{"m", cdim_pckd},
				&la.IOpt{"n", p}, la.OptTrans)
			lapack.Trtrs(QA, y, la.OptUpper, &la.IOpt{"n", p})
			//checkpnt.Check("80solve_qr", minor)

			unpack(u, z, dims)
			checkpnt.Check("90solve_qr", minor)
			return nil
		}
		return solve, err
	}
	return factor, nil
}

//    Solution of KKT equations by reduction to a 2 x 2 system, a QR
//    factorization to eliminate the equality constraints, and a dense
//    Cholesky factorization of order n-p.
//
//    Computes the QR factorization
//
//        A' = [Q1, Q2] * [R; 0]
//
//    and returns a function that (1) computes the Cholesky factorization
//
//        Q_2^T * (H + GG^T * W^{-1} * W^{-T} * GG) * Q2 = L * L^T,
//
//    given H, Df, W, where GG = [Df; G], and (2) returns a function for
//    solving
//
//        [ H    A'   GG'    ]   [ ux ]   [ bx ]
//        [ A    0    0      ] * [ uy ] = [ by ].
//        [ GG   0    -W'*W  ]   [ uz ]   [ bz ]
//
//    H is n x n,  A is p x n, Df is mnl x n, G is N x n where
//    N = dims['l'] + sum(dims['q']) + sum( k**2 for k in dims['s'] ).
//
func kktChol(G *matrix.FloatMatrix, dims *sets.DimensionSet, A *matrix.FloatMatrix, mnl int) (kktFactor, error) {

	p, n := A.Size()
	cdim := mnl + dims.Sum("l", "q") + dims.SumSquared("s")
	cdim_pckd := mnl + dims.Sum("l", "q") + dims.SumPacked("s")

	QA := A.Transpose()
	tauA := matrix.FloatZeros(p, 1)
	lapack.Geqrf(QA, tauA)

	Gs := matrix.FloatZeros(cdim, n)
	K := matrix.FloatZeros(n, n)
	bzp := matrix.FloatZeros(cdim_pckd, 1)
	yy := matrix.FloatZeros(p, 1)
	checkpnt.AddMatrixVar("tauA", tauA)
	checkpnt.AddMatrixVar("Gs", Gs)
	checkpnt.AddMatrixVar("K", K)

	factor := func(W *sets.FloatMatrixSet, H, Df *matrix.FloatMatrix) (KKTFunc, error) {
		// Compute
		//
		//     K = [Q1, Q2]' * (H + GG' * W^{-1} * W^{-T} * GG) * [Q1, Q2]
		//
		// and take the Cholesky factorization of the 2,2 block
		//
		//     Q_2' * (H + GG^T * W^{-1} * W^{-T} * GG) * Q2.

		var err error = nil
		minor := 0
		if !checkpnt.MinorEmpty() {
			minor = checkpnt.MinorTop()
		}
		// Gs = W^{-T} * GG in packed storage.
		if mnl > 0 {
			Gs.SetSubMatrix(0, 0, Df)
		}
		Gs.SetSubMatrix(mnl, 0, G)
		checkpnt.Check("00factor_chol", minor)
		scale(Gs, W, true, true)
		pack2(Gs, dims, mnl)
		//checkpnt.Check("10factor_chol", minor)

		// K = [Q1, Q2]' * (H + Gs' * Gs) * [Q1, Q2].
		blas.SyrkFloat(Gs, K, 1.0, 0.0, la.OptTrans, &la.IOpt{"k", cdim_pckd})
		if H != nil {
			K.SetSubMatrix(0, 0, matrix.Plus(H, K.GetSubMatrix(0, 0, H.Rows(), H.Cols())))
		}
		//checkpnt.Check("20factor_chol", minor)
		symm(K, n, 0)
		lapack.Ormqr(QA, tauA, K, la.OptLeft, la.OptTrans)
		lapack.Ormqr(QA, tauA, K, la.OptRight)
		//checkpnt.Check("30factor_chol", minor)

		// Cholesky factorization of 2,2 block of K.
		lapack.Potrf(K, &la.IOpt{"n", n - p}, &la.IOpt{"offseta", p * (n + 1)})
		checkpnt.Check("40factor_chol", minor)

		solve := func(x, y, z *matrix.FloatMatrix) (err error) {
			// Solve
			//
			//     [ 0          A'  GG'*W^{-1} ]   [ ux   ]   [ bx        ]
			//     [ A          0   0          ] * [ uy   ] = [ by        ]
			//     [ W^{-T}*GG  0   -I         ]   [ W*uz ]   [ W^{-T}*bz ]
			//
			// and return ux, uy, W*uz.
			//
			// On entry, x, y, z contain bx, by, bz.  On exit, they contain
			// the solution ux, uy, W*uz.
			//
			// If we change variables ux = Q1*v + Q2*w, the system becomes
			//
			//     [ K11 K12 R ]   [ v  ]   [Q1'*(bx+GG'*W^{-1}*W^{-T}*bz)]
			//     [ K21 K22 0 ] * [ w  ] = [Q2'*(bx+GG'*W^{-1}*W^{-T}*bz)]
			//     [ R^T 0   0 ]   [ uy ]   [by                           ]
			//
			//     W*uz = W^{-T} * ( GG*ux - bz ).
			minor := 0
			if !checkpnt.MinorEmpty() {
				minor = checkpnt.MinorTop()
			}

			// bzp := W^{-T} * bz in packed storage
			scale(z, W, true, true)
			pack(z, bzp, dims, &la.IOpt{"mnl", mnl})

			// x := [Q1, Q2]' * (x + Gs' * bzp)
			//    = [Q1, Q2]' * (bx + Gs' * W^{-T} * bz)
			blas.GemvFloat(Gs, bzp, x, 1.0, 1.0, la.OptTrans, &la.IOpt{"m", cdim_pckd})
			lapack.Ormqr(QA, tauA, x, la.OptLeft, la.OptTrans)

			// y := x[:p]
			//    = Q1' * (bx + Gs' * W^{-T} * bz)
			blas.Copy(y, yy)
			blas.Copy(x, y, &la.IOpt{"n", p})

			// x[:p] := v = R^{-T} * by
			blas.Copy(yy, x)
			lapack.Trtrs(QA, x, la.OptUpper, la.OptTrans, &la.IOpt{"n", p})

			// x[p:] := K22^{-1} * (x[p:] - K21*x[:p])
			//        = K22^{-1} * (Q2' * (bx + Gs' * W^{-T} * bz) - K21*v)
			blas.GemvFloat(K, x, x, -1.0, 1.0, &la.IOpt{"m", n - p}, &la.IOpt{"n", p},
				&la.IOpt{"offseta", p}, &la.IOpt{"offsety", p})
			lapack.Potrs(K, x, &la.IOpt{"n", n - p}, &la.IOpt{"offseta", p * (n + 1)},
				&la.IOpt{"offsetb", p})

			// y := y - [K11, K12] * x
			//    = Q1' * (bx + Gs' * W^{-T} * bz) - K11*v - K12*w
			blas.GemvFloat(K, x, y, -1.0, 1.0, &la.IOpt{"m", p}, &la.IOpt{"n", n})

			// y := R^{-1}*y
			//    = R^{-1} * (Q1' * (bx + Gs' * W^{-T} * bz) - K11*v
			//      - K12*w)
			lapack.Trtrs(QA, y, la.OptUpper, &la.IOpt{"n", p})

			// x := [Q1, Q2] * x
			lapack.Ormqr(QA, tauA, x, la.OptLeft)

			// bzp := Gs * x - bzp.
			//      = W^{-T} * ( GG*ux - bz ) in packed storage.
			// Unpack and copy to z.
			blas.GemvFloat(Gs, x, bzp, 1.0, -1.0, &la.IOpt{"m", cdim_pckd})
			unpack(bzp, z, dims, &la.IOpt{"mnl", mnl})

			checkpnt.Check("90solve_chol", minor)
			return nil
		}
		return solve, err
	}
	return factor, nil
}

/*
   Solution of KKT equations by reduction to a 2 x 2 system, a sparse
   or dense Cholesky factorization of order n to eliminate the 1,1
   block, and a sparse or dense Cholesky factorization of order p.
   Implemented only for problems with no second-order or semidefinite
   cone constraints.

   Returns a function that (1) computes Cholesky factorizations of
   the matrices

       S = H + GG' * W^{-1} * W^{-T} * GG,
       K = A * S^{-1} *A'

   or (if K is singular in the first call to the function), the matrices

       S = H + GG' * W^{-1} * W^{-T} * GG + A' * A,
       K = A * S^{-1} * A',

   given H, Df, W, where GG = [Df; G], and (2) returns a function for
   solving

       [ H     A'   GG'   ]   [ ux ]   [ bx ]
       [ A     0    0     ] * [ uy ] = [ by ].
       [ GG    0   -W'*W  ]   [ uz ]   [ bz ]

   H is n x n,  A is p x n, Df is mnl x n, G is dims['l'] x n.

*/

type chol2Data struct {
	firstcall     bool
	singular      bool
	A, G          *matrix.FloatMatrix
	dims          *sets.DimensionSet
	Gs, S, K, Dfs *matrix.FloatMatrix
}

func kktChol2(G *matrix.FloatMatrix, dims *sets.DimensionSet, A *matrix.FloatMatrix, mnl int) (kktFactor, error) {

	if len(dims.At("q")) > 0 || len(dims.At("s")) > 0 {
		return nil, errors.New("'chol2' solver only for problems with no second-order or " +
			"semidefinite cone constraints")
	}

	p, n := A.Size()
	ml := dims.At("l")[0]
	F := &chol2Data{firstcall: true, singular: false, A: A, G: G, dims: dims}

	factor := func(W *sets.FloatMatrixSet, H, Df *matrix.FloatMatrix) (KKTFunc, error) {
		var err error = nil
		minor := 0
		if !checkpnt.MinorEmpty() {
			minor = checkpnt.MinorTop()
		}
		if F.firstcall {
			F.Gs = matrix.FloatZeros(F.G.Size())
			if mnl > 0 {
				F.Dfs = matrix.FloatZeros(Df.Size())
			}
			F.S = matrix.FloatZeros(n, n)
			F.K = matrix.FloatZeros(p, p)
			checkpnt.AddMatrixVar("Gs", F.Gs)
			checkpnt.AddMatrixVar("Dfs", F.Dfs)
			checkpnt.AddMatrixVar("S", F.S)
			checkpnt.AddMatrixVar("K", F.K)
		}

		if mnl > 0 {
			dnli := matrix.FloatZeros(mnl, mnl)
			dnli.SetIndexesFromArray(W.At("dnli")[0].FloatArray(), matrix.DiagonalIndexes(dnli)...)
			blas.GemmFloat(dnli, Df, F.Dfs, 1.0, 0.0)
		}
		checkpnt.Check("02factor_chol2", minor)
		di := matrix.FloatZeros(ml, ml)
		di.SetIndexesFromArray(W.At("di")[0].FloatArray(), matrix.DiagonalIndexes(di)...)
		err = blas.GemmFloat(di, G, F.Gs, 1.0, 0.0)
		checkpnt.Check("06factor_chol2", minor)

		if F.firstcall {
			blas.SyrkFloat(F.Gs, F.S, 1.0, 0.0, la.OptTrans)
			if mnl > 0 {
				blas.SyrkFloat(F.Dfs, F.S, 1.0, 1.0, la.OptTrans)
			}
			if H != nil {
				F.S.Plus(H)
			}
			checkpnt.Check("10factor_chol2", minor)
			err = lapack.Potrf(F.S)
			if err != nil {
				err = nil // reset error
				F.singular = true
				// original code recreates F.S as dense if it is sparse and
				// A is dense, we don't do it as currently no sparse matrices
				//F.S = matrix.FloatZeros(n, n)
				//checkpnt.AddMatrixVar("S", F.S)
				blas.SyrkFloat(F.Gs, F.S, 1.0, 0.0, la.OptTrans)
				if mnl > 0 {
					blas.SyrkFloat(F.Dfs, F.S, 1.0, 1.0, la.OptTrans)
				}
				checkpnt.Check("14factor_chol2", minor)
				blas.SyrkFloat(F.A, F.S, 1.0, 1.0, la.OptTrans)
				if H != nil {
					F.S.Plus(H)
				}
				lapack.Potrf(F.S)
			}
			F.firstcall = false
			checkpnt.Check("20factor_chol2", minor)
		} else {
			blas.SyrkFloat(F.Gs, F.S, 1.0, 0.0, la.OptTrans)
			if mnl > 0 {
				blas.SyrkFloat(F.Dfs, F.S, 1.0, 1.0, la.OptTrans)
			}
			if H != nil {
				F.S.Plus(H)
			}
			checkpnt.Check("40factor_chol2", minor)
			if F.singular {
				blas.SyrkFloat(F.A, F.S, 1.0, 1.0, la.OptTrans)
			}
			lapack.Potrf(F.S)
			checkpnt.Check("50factor_chol2", minor)
		}

		// Asct := L^{-1}*A'.  Factor K = Asct'*Asct.
		Asct := F.A.Transpose()
		blas.TrsmFloat(F.S, Asct, 1.0)
		blas.SyrkFloat(Asct, F.K, 1.0, 0.0, la.OptTrans)
		lapack.Potrf(F.K)
		checkpnt.Check("90factor_chol2", minor)

		solve := func(x, y, z *matrix.FloatMatrix) (err error) {
			// Solve
			//
			//     [ H          A'  GG'*W^{-1} ]   [ ux   ]   [ bx        ]
			//     [ A          0   0          ] * [ uy   ] = [ by        ]
			//     [ W^{-T}*GG  0   -I         ]   [ W*uz ]   [ W^{-T}*bz ]
			//
			// and return ux, uy, W*uz.
			//
			// If not F['singular']:
			//
			//     K*uy = A * S^{-1} * ( bx + GG'*W^{-1}*W^{-T}*bz ) - by
			//     S*ux = bx + GG'*W^{-1}*W^{-T}*bz - A'*uy
			//     W*uz = W^{-T} * ( GG*ux - bz ).
			//
			// If F['singular']:
			//
			//     K*uy = A * S^{-1} * ( bx + GG'*W^{-1}*W^{-T}*bz + A'*by )
			//            - by
			//     S*ux = bx + GG'*W^{-1}*W^{-T}*bz + A'*by - A'*y.
			//     W*uz = W^{-T} * ( GG*ux - bz ).

			minor := 0
			if !checkpnt.MinorEmpty() {
				minor = checkpnt.MinorTop()
			}

			// z := W^{-1} * z = W^{-1} * bz
			scale(z, W, true, true)
			checkpnt.Check("10solve_chol2", minor)

			// If not F['singular']:
			//     x := L^{-1} * P * (x + GGs'*z)
			//        = L^{-1} * P * (x + GG'*W^{-1}*W^{-T}*bz)
			//
			// If F['singular']:
			//     x := L^{-1} * P * (x + GGs'*z + A'*y))
			//        = L^{-1} * P * (x + GG'*W^{-1}*W^{-T}*bz + A'*y)
			if mnl > 0 {
				blas.GemvFloat(F.Dfs, z, x, 1.0, 1.0, la.OptTrans)
			}
			blas.GemvFloat(F.Gs, z, x, 1.0, 1.0, la.OptTrans, &la.IOpt{"offsetx", mnl})
			//checkpnt.Check("20solve_chol2", minor)
			if F.singular {
				blas.GemvFloat(F.A, y, x, 1.0, 1.0, la.OptTrans)
			}
			checkpnt.Check("30solve_chol2", minor)
			blas.TrsvFloat(F.S, x)
			//checkpnt.Check("50solve_chol2", minor)

			// y := K^{-1} * (Asc*x - y)
			//    = K^{-1} * (A * S^{-1} * (bx + GG'*W^{-1}*W^{-T}*bz) - by)
			//      (if not F['singular'])
			//    = K^{-1} * (A * S^{-1} * (bx + GG'*W^{-1}*W^{-T}*bz +
			//      A'*by) - by)
			//      (if F['singular']).
			blas.GemvFloat(Asct, x, y, 1.0, -1.0, la.OptTrans)
			//checkpnt.Check("55solve_chol2", minor)
			lapack.Potrs(F.K, y)
			//checkpnt.Check("60solve_chol2", minor)

			// x := P' * L^{-T} * (x - Asc'*y)
			//    = S^{-1} * (bx + GG'*W^{-1}*W^{-T}*bz - A'*y)
			//      (if not F['singular'])
			//    = S^{-1} * (bx + GG'*W^{-1}*W^{-T}*bz + A'*by - A'*y)
			//      (if F['singular'])
			blas.GemvFloat(Asct, y, x, -1.0, 1.0)
			blas.TrsvFloat(F.S, x, la.OptTrans)
			checkpnt.Check("70solve_chol2", minor)

			// W*z := GGs*x - z = W^{-T} * (GG*x - bz)
			if mnl > 0 {
				blas.GemvFloat(F.Dfs, x, z, 1.0, -1.0)
			}
			blas.GemvFloat(F.Gs, x, z, 1.0, -1.0, &la.IOpt{"offsety", mnl})

			checkpnt.Check("90solve_chol2", minor)
			return nil
		}
		return solve, err
	}
	return factor, nil
}

// Local Variables:
// tab-width: 4
// indent-tabs-mode: nil
// End:
