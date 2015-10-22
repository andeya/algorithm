// Copyright (c) Harri Rautila, 2012,2013

// This file is part of github.com/hrautila/matops package. It is free software,
// distributed under the terms of GNU Lesser General Public License Version 3, or
// any later version. See the COPYING tile included in this archive.

package matops

import (
	"errors"
	"github.com/henrylee2cn/algorithm/matops/calgo"
	"github.com/henrylee2cn/algorithm/matrix"
	"math"
	//"fmt"
)

// Compute
//      Y = alpha*A*X + beta*Y
//      Y = alpha*A.T*X + beta*Y  ; flags = TRANSA
//
//    A is M*N or N*M generic matrix,
//    X is row or column vector of length N
//    Y is row or column vector of legth M.
//
// MVMult is vector orientation agnostic. It does not matter if Y, X are row or
// column vectors, they are always handled as if they were column vectors.
func MVMult(Y, A, X *matrix.FloatMatrix, alpha, beta float64, flags Flags) error {

	if A.Rows() == 0 || A.Cols() == 0 {
		return nil
	}
	if Y.Rows() != 1 && Y.Cols() != 1 {
		return errors.New("Y not a vector.")
	}
	if X.Rows() != 1 && X.Cols() != 1 {
		return errors.New("X not a vector.")
	}

	Ar := A.FloatArray()
	ldA := A.LeadingIndex()
	Yr := Y.FloatArray()
	incY := 1
	lenY := Y.NumElements()
	if Y.Rows() == 1 {
		// row vector
		incY = Y.LeadingIndex()
	}
	Xr := X.FloatArray()
	incX := 1
	lenX := X.NumElements()
	if X.Rows() == 1 {
		// row vector
		incX = X.LeadingIndex()
	}
	// NOTE: This could diveded to parallel tasks by rows.
	calgo.DMultMV(Yr, Ar, Xr, alpha, beta, calgo.Flags(flags), incY, ldA, incX,
		0, lenX, 0, lenY, vpLen, mB)
	return nil
}

// Matrix-vector rank update A = A + alpha*X*Y.T
//    A is M*N generic matrix,
//    X is row or column vector of length M
//    Y is row or column vector of legth N.
func MVRankUpdate(A, X, Y *matrix.FloatMatrix, alpha float64) error {

	if A.Rows() == 0 || A.Cols() == 0 {
		return nil
	}
	if Y.Rows() != 1 && Y.Cols() != 1 {
		return errors.New("Y not a vector.")
	}
	if X.Rows() != 1 && X.Cols() != 1 {
		return errors.New("X not a vector.")
	}

	Ar := A.FloatArray()
	ldA := A.LeadingIndex()
	Yr := Y.FloatArray()
	incY := 1
	if Y.Rows() == 1 {
		// row vector
		incY = Y.LeadingIndex()
	}
	Xr := X.FloatArray()
	incX := 1
	if X.Rows() == 1 {
		// row vector
		incX = X.LeadingIndex()
	}
	// NOTE: This could diveded to parallel tasks like matrix-matrix multiplication
	calgo.DRankMV(Ar, Xr, Yr, alpha, ldA, incX, incY, 0, A.Cols(), 0, A.Rows(), 0, 0)
	return nil
}

// Matrix-vector symmetric rank update A = A + alpha*X*X.T
//   A is N*N symmetric,
//   X is row or column vector of length N.
func MVRankUpdateSym(A, X *matrix.FloatMatrix, alpha float64, flags Flags) error {

	if A.Rows() == 0 || A.Cols() == 0 {
		return nil
	}
	if X.Rows() != 1 && X.Cols() != 1 {
		return errors.New("X not a vector.")
	}

	Ar := A.FloatArray()
	ldA := A.LeadingIndex()
	Xr := X.FloatArray()
	incX := 1
	if X.Rows() == 1 {
		// row vector
		incX = X.LeadingIndex()
	}
	// NOTE: This could diveded to parallel tasks per column.
	calgo.DSymmRankMV(Ar, Xr, alpha, calgo.Flags(flags), ldA, incX, 0, A.Cols(), 0)
	return nil
}

// Matrix-vector symmetric rank 2 update A = A + alpha*X*Y.T + alpha*X.T*Y
//   A is N*N symmetric matrix,
//   X is row or column vector of length N
//   Y is row or column vector of legth N.
func MVRankUpdate2Sym(A, X, Y *matrix.FloatMatrix, alpha float64, flags Flags) error {

	if A.Rows() == 0 || A.Cols() == 0 {
		return nil
	}
	if Y.Rows() != 1 && Y.Cols() != 1 {
		return errors.New("Y not a vector.")
	}
	if X.Rows() != 1 && X.Cols() != 1 {
		return errors.New("X not a vector.")
	}

	Ar := A.FloatArray()
	ldA := A.LeadingIndex()
	Yr := Y.FloatArray()
	incY := 1
	if Y.Rows() == 1 {
		// row vector
		incY = Y.LeadingIndex()
	}
	Xr := X.FloatArray()
	incX := 1
	if X.Rows() == 1 {
		// row vector
		incX = X.LeadingIndex()
	}
	// NOTE: This could diveded to parallel tasks like matrix-matrix multiplication
	calgo.DSymmRank2MV(Ar, Xr, Yr, alpha, calgo.Flags(flags), ldA, incY, incX, 0, A.Cols(), 0)
	return nil
}

// Matrix-vector triangular update A = A + alpha*X*Y.T
//   A is N*N matrix,
//   X is row or column vector of length N
//   Y is row or column vector of legth N.
//   flags is UPPER or LOWER
func MVUpdateTrm(A, X, Y *matrix.FloatMatrix, alpha float64, flags Flags) error {

	if A.Rows() == 0 || A.Cols() == 0 {
		return nil
	}
	if Y.Rows() != 1 && Y.Cols() != 1 {
		return errors.New("Y not a vector.")
	}
	if X.Rows() != 1 && X.Cols() != 1 {
		return errors.New("X not a vector.")
	}

	Ar := A.FloatArray()
	ldA := A.LeadingIndex()
	Yr := Y.FloatArray()
	incY := 1
	if Y.Rows() == 1 {
		// row vector
		incY = Y.LeadingIndex()
	}
	Xr := X.FloatArray()
	incX := 1
	if X.Rows() == 1 {
		// row vector
		incX = X.LeadingIndex()
	}
	// NOTE: This could diveded to parallel tasks like matrix-matrix multiplication
	calgo.DTrmUpdMV(Ar, Xr, Yr, alpha, calgo.Flags(flags), ldA, incX, incY, 0, A.Cols(), nB)
	return nil
}

// Matrix-vector solve X = A.-1*X or X = A.-T*X
//   A is N*N tridiagonal lower or upper,
//   X is row or column vector of length N.
// flags
//   LOWER  A is lower tridiagonal
//   UPPER  A is upper tridiagonal
//   UNIT   A diagonal is unit
//   TRANSA A is transpose
func MVSolveTrm(X, A *matrix.FloatMatrix, alpha float64, flags Flags) error {

	if A.Rows() == 0 || A.Cols() == 0 {
		return nil
	}
	if X.Rows() != 1 && X.Cols() != 1 {
		return errors.New("X not a vector.")
	}

	Ar := A.FloatArray()
	ldA := A.LeadingIndex()
	Xr := X.FloatArray()
	incX := 1
	if X.Rows() == 1 {
		// row vector
		incX = X.LeadingIndex()
	}
	calgo.DSolveBlkMV(Xr, Ar, calgo.Flags(flags), incX, ldA, A.Cols(), nB)
	return nil
}

// Tridiagonal multiplication; X = A*X
func MVMultTrm(X, A *matrix.FloatMatrix, flags Flags) error {

	if A.Rows() == 0 || A.Cols() == 0 {
		return nil
	}
	if X.Rows() != 1 && X.Cols() != 1 {
		return errors.New("X not a vector.")
	}

	Ar := A.FloatArray()
	ldA := A.LeadingIndex()
	Xr := X.FloatArray()
	incX := 1
	if X.Rows() == 1 {
		// row vector
		incX = X.LeadingIndex()
	}
	calgo.DTrimvUnblkMV(Xr, Ar, calgo.Flags(flags), incX, ldA, A.Cols())
	return nil
}

// DiffNorm2: sqrt(||X - Y||^2)
func DiffNorm2(X, Y *matrix.FloatMatrix) float64 {
	if X == nil || Y == nil {
		return math.NaN()
	}
	if !isVector(X) {
		return math.NaN()
	}
	if !isVector(Y) {
		return math.NaN()
	}
	Xr := X.FloatArray()
	incX := 1
	if X.Cols() != 1 {
		// Row vector
		incX = X.LeadingIndex()
	}
	Yr := Y.FloatArray()
	incY := 1
	if Y.Cols() != 1 {
		// Row vector
		incY = Y.LeadingIndex()
	}
	return calgo.DiffNorm2(Xr, Yr, incX, incY, X.NumElements())
}

// Inner product: alpha * X * Y
func Dot(X, Y *matrix.FloatMatrix, alpha float64) float64 {
	if X == nil || Y == nil {
		return math.NaN()
	}
	if !isVector(X) {
		return math.NaN()
	}
	if !isVector(Y) {
		return math.NaN()
	}
	Xr := X.FloatArray()
	incX := 1
	if X.Cols() != 1 {
		// Row vector
		incX = X.LeadingIndex()
	}
	Yr := Y.FloatArray()
	incY := 1
	if Y.Cols() != 1 {
		// Row vector
		incY = Y.LeadingIndex()
	}
	return calgo.DDot(Xr, Yr, alpha, incX, incY, X.NumElements())
}

// Add inner product: Z[index] = beta*Z[index] + alpha * X * Y
func AddDot(Z, X, Y *matrix.FloatMatrix, alpha, beta float64, index int) {
	if X == nil || Y == nil {
		return
	}
	if X.NumElements() == 0 || Y.NumElements() == 0 {
		return
	}
	if !isVector(X) {
		return
	}
	if !isVector(Y) {
		return
	}
	Xr := X.FloatArray()
	incX := 1
	if X.Cols() != 1 {
		// Row vector
		incX = X.LeadingIndex()
	}
	Yr := Y.FloatArray()
	incY := 1
	if Y.Cols() != 1 {
		// Row vector
		incY = Y.LeadingIndex()
	}
	Zr := Z.FloatArray()
	incZ := 1
	if Z.Cols() != 1 {
		// Row vector
		incZ = Z.LeadingIndex()
	}
	calgo.DDotSum(Zr[incZ*index:], Xr, Yr, alpha, beta, incZ, incX, incY, X.NumElements())
}

// Y := alpha * X + Y
func Axpy(Y, X *matrix.FloatMatrix, alpha float64) {
	if X == nil || Y == nil {
		return
	}
	if !isVector(X) {
		return
	}
	if !isVector(Y) {
		return
	}
	Xr := X.FloatArray()
	incX := 1
	if X.Cols() != 1 {
		// Row vector
		incX = X.LeadingIndex()
	}
	Yr := Y.FloatArray()
	incY := 1
	if Y.Cols() != 1 {
		// Row vector
		incY = Y.LeadingIndex()
	}
	calgo.DAxpy(Xr, Yr, alpha, incX, incY, X.NumElements())
	return
}

// Norm2 of vector: sqrt(||x||^2)
func Norm2(X *matrix.FloatMatrix) float64 {
	if X == nil {
		return math.NaN()
	}
	if !isVector(X) {
		return math.NaN()
	}
	Xr := X.FloatArray()
	incX := 1
	if X.Cols() != 1 {
		// Row vector
		incX = X.LeadingIndex()
	}
	return calgo.DNorm2(Xr, incX, X.NumElements())
}

// sum(|x|)
func ASum(X *matrix.FloatMatrix) float64 {
	if X == nil {
		return math.NaN()
	}
	if !isVector(X) {
		return math.NaN()
	}
	Xr := X.FloatArray()
	incX := 1
	if X.Cols() != 1 {
		// Row vector
		incX = X.LeadingIndex()
	}
	return calgo.DAsum(Xr, incX, X.NumElements())
}

// max |x|
func AMax(X *matrix.FloatMatrix) float64 {
	ix := IAMax(X)
	if ix == -1 {
		return math.NaN()
	}
	return X.GetIndex(ix)
}

// index of max |x|
func IAMax(X *matrix.FloatMatrix) int {
	if X == nil {
		return -1
	}
	if !isVector(X) {
		return -1
	}
	Xr := X.FloatArray()
	incX := 1
	if X.Cols() != 1 {
		// Row vector
		incX = X.LeadingIndex()
	}
	return calgo.DIAMax(Xr, incX, X.NumElements())
}

// Inverse scaling of vector. X = X / alpha.
func InvScale(X *matrix.FloatMatrix, alpha float64) {
	if X == nil || X.NumElements() == 0 {
		return
	}
	if !isVector(X) {
		return
	}
	Xr := X.FloatArray()
	incX := 1
	if X.Cols() != 1 {
		// Row vector
		incX = X.LeadingIndex()
	}
	calgo.DInvScal(Xr, alpha, incX, X.NumElements())
}

// Scaling with scalar: X = alpha * X
func Scale(X *matrix.FloatMatrix, alpha float64) {
	if X == nil || X.NumElements() == 0 {
		return
	}
	if !isVector(X) {
		return
	}
	Xr := X.FloatArray()
	incX := 1
	if X.Cols() != 1 {
		// Row vector
		incX = X.LeadingIndex()
	}
	calgo.DScal(Xr, alpha, incX, X.NumElements())
}

// Swap X and Y.
func Swap(X, Y *matrix.FloatMatrix) {
	if X == nil || Y == nil {
		return
	}
	if X.NumElements() == 0 || Y.NumElements() == 0 {
		return
	}
	if !isVector(X) {
		return
	}
	if !isVector(Y) {
		return
	}
	Xr := X.FloatArray()
	incX := 1
	if X.Cols() != 1 {
		// Row vector
		incX = X.LeadingIndex()
	}
	Yr := Y.FloatArray()
	incY := 1
	if Y.Cols() != 1 {
		// Row vector
		incY = Y.LeadingIndex()
	}
	calgo.DSwap(Xr, Yr, incX, incY, X.NumElements())
}

// Local Variables:
// tab-width: 4
// indent-tabs-mode: nil
// End:
