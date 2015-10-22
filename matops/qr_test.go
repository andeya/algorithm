// Copyright (c) Harri Rautila, 2012,2013

// This file is part of github.com/hrautila/matops package. It is free software,
// distributed under the terms of GNU Lesser General Public License Version 3, or
// any later version. See the COPYING tile included in this archive.

package matops

import (
	"github.com/henrylee2cn/algorithm/matrix"
	"testing"
)

func TestDecomposeQR(t *testing.T) {
	M := 60
	N := 40
	nb := 4
	A := matrix.FloatUniform(M, N)
	W := matrix.FloatZeros(N, nb)
	tau := matrix.FloatZeros(N, 1)

	// QR = Q*R
	QR, _ := DecomposeQR(A.Copy(), tau, W, nb)

	// A2 = Q*R
	A2 := TriU(QR.Copy())
	MultQ(A2, QR, tau, W, LEFT, 0)

	A.Minus(A2)
	// ||A - Q*R||_1
	nrm := NormP(A, NORM_ONE)
	t.Logf("||A - Q*R||_1: %e\n", nrm)
}

func TestDecomposeQRT(t *testing.T) {
	M := 60
	N := 40
	nb := 12
	A := matrix.FloatUniform(M, N)
	W := matrix.FloatZeros(N, N)
	T := matrix.FloatZeros(N, N)

	// QR = Q*R
	QR, _ := DecomposeQRT(A.Copy(), T, W, nb)
	// A2 = Q*R
	A2 := TriU(QR.Copy())
	MultQT(A2, QR, T, W, LEFT, 0)

	A.Minus(A2)
	// ||A - Q*R||_1
	nrm := NormP(A, NORM_ONE)
	t.Logf("||A - Q*R||_1: %e\n", nrm)
}

func TestMultQ(t *testing.T) {
	M := 60
	N := 40
	K := 30
	nb := 12
	A := matrix.FloatUniform(M, N)
	B := matrix.FloatUniform(M, K)
	W := matrix.FloatZeros(N, nb)
	tau := matrix.FloatZeros(N, 1)
	X := B.Copy()

	// QR = Q*R
	QR, err := DecomposeQR(A.Copy(), tau, W, nb)
	if err != nil {
		t.Logf("decompose-err: %v\n", err)
	}

	// B - Q*Q.T*B = 0
	MultQ(X, QR, tau, W, LEFT|TRANS, nb)
	MultQ(X, QR, tau, W, LEFT, nb)
	// X = Q*Q.T*B

	B.Minus(X)
	// ||B - Q*Q.T*B||_1
	nrm := NormP(B, NORM_ONE)
	t.Logf("||B - Q*Q.T*B||_1: %e\n", nrm)
}

func TestMultQT(t *testing.T) {
	M := 60
	N := 40
	K := 30
	nb := 12
	A := matrix.FloatUniform(M, N)
	B := matrix.FloatUniform(M, K)
	W := matrix.FloatZeros(N, nb)
	T := matrix.FloatZeros(N, N)

	// QR = Q*R
	QR, err := DecomposeQRT(A.Copy(), T, W, nb)
	if err != nil {
		t.Logf("decompose-err: %v\n", err)
	}

	// compute: B - Q*Q.T*B = 0

	// X = Q*Q.T*B
	X := B.Copy()
	MultQT(X, QR, T, W, LEFT|TRANS, nb)
	MultQT(X, QR, T, W, LEFT, nb)
	B.Minus(X)

	// ||B - Q*Q.T*B||_1
	nrm := NormP(B, NORM_ONE)
	t.Logf("||B - Q*Q.T*B||_1: %e\n", nrm)
}

func TestSolveLeastSquaresQR(t *testing.T) {
	M := 60
	N := 40
	K := 30
	nb := 12

	A := matrix.FloatUniform(M, N)
	B := matrix.FloatZeros(M, K)
	X0 := matrix.FloatUniform(N, K)

	// B = A*X0
	Mult(B, A, X0, 1.0, 0.0, NOTRANS)
	W := matrix.FloatZeros(N, nb)
	tau := matrix.FloatZeros(N, 1)

	QR, err := DecomposeQR(A, tau, W, nb)
	if err != nil {
		t.Logf("decompose error: %v\n", err)
	}
	// B' = A.-1*B
	err = SolveQR(B, QR, tau, W, NOTRANS, nb)

	// expect B[0:N, 0:K] == X0, B[N:M, 0:K] == 0.0
	var Xref matrix.FloatMatrix
	Bref := matrix.FloatZeros(M, K)
	Bref.SubMatrix(&Xref, 0, 0, N, K)
	Xref.Plus(X0)
	Bref.Minus(B)

	t.Logf("||B - A*X0||_1: %e\n", NormP(Bref, NORM_ONE))
}

func TestSolveLeastSquaresQRT(t *testing.T) {
	M := 60
	N := 40
	K := 30
	nb := 12

	A := matrix.FloatUniform(M, N)
	B := matrix.FloatZeros(M, K)
	X0 := matrix.FloatUniform(N, K)

	// B = A*X0
	Mult(B, A, X0, 1.0, 0.0, NOTRANS)
	W := matrix.FloatZeros(N, nb)
	T := matrix.FloatZeros(N, N)

	QR, err := DecomposeQRT(A, T, W, nb)
	if err != nil {
		t.Logf("decompose error: %v\n", err)
	}
	// B' = A.-1*B
	err = SolveQRT(B, QR, T, W, NOTRANS, nb)

	// expect B[0:N, 0:K] == X0, B[N:M, 0:K] == 0.0
	var Xref matrix.FloatMatrix
	Bref := matrix.FloatZeros(M, K)
	Bref.SubMatrix(&Xref, 0, 0, N, K)
	Xref.Plus(X0)
	Bref.Minus(B)
	t.Logf("\nmin ||B - A*X0||\n\twhere B = A*X0\n")
	t.Logf("||B - A*X0||_1 ~~ 0.0: %e\n", NormP(Bref, NORM_ONE))
}

func TestSolveMinQRT(t *testing.T) {
	var QR *matrix.FloatMatrix
	var err error
	N := 40
	K := 30
	nb := 12

	A := matrix.FloatUniform(N, N)
	B0 := matrix.FloatUniform(N, K)
	W := matrix.FloatZeros(N, nb)
	T := matrix.FloatZeros(N, N)

	QR, err = DecomposeQRT(A.Copy(), T, W, nb)
	if err != nil {
		t.Logf("decompose error: %v\n", err)
	}
	// B' = A.-1*B
	B := B0.Copy()
	err = SolveQRT(B, QR, T, W, TRANS, nb)
	Mult(B0, A, B, -1.0, 1.0, TRANSA)
	t.Logf("\nmin ||X||\n\ts.t. A.T*X = B\n\twhere A is N-by-N, B is N-by-K\n")
	t.Logf("||B - A.T*X||_1 ~~ 0.0: %e\n", NormP(B0, NORM_ONE))
}

// Local Variables:
// tab-width: 4
// indent-tabs-mode: nil
// End:
