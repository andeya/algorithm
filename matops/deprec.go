package matops

import (
	"errors"
	"github.com/henrylee2cn/algorithm/matops/calgo"
	"github.com/henrylee2cn/algorithm/matrix"
)

func Mult0(C, A, B *matrix.FloatMatrix, alpha, beta float64, flags Flags) error {
	if A.Cols() != B.Rows() {
		return errors.New("A.cols != B.rows: size mismatch")
	}
	psize := int64(C.NumElements()) * int64(A.Cols())
	Ar := A.FloatArray()
	ldA := A.LeadingIndex()
	Br := B.FloatArray()
	ldB := B.LeadingIndex()
	Cr := C.FloatArray()
	ldC := C.LeadingIndex()

	if nWorker <= 1 || psize <= limitOne {
		calgo.DMult0(Cr, Ar, Br, alpha, beta, calgo.Flags(flags), ldC, ldA, ldB, B.Rows(),
			0, C.Cols(), 0, C.Rows(),
			vpLen, nB, mB)
		return nil
	}
	// here we have more than one worker available
	worker := func(cstart, cend, rstart, rend int, ready chan int) {
		calgo.DMult0(Cr, Ar, Br, alpha, beta, calgo.Flags(flags), ldC, ldA, ldB, B.Rows(),
			cstart, cend, rstart, rend, vpLen, nB, mB)
		ready <- 1
	}
	colworks, rowworks := divideWork(C.Rows(), C.Cols(), nWorker)
	scheduleWork(colworks, rowworks, C.Cols(), C.Rows(), worker)
	return nil
}

// Calculate C = alpha*A*B + beta*C, C is M*N, A is M*P and B is P*N
func MMMultNoTrans(C, A, B *matrix.FloatMatrix, alpha, beta float64) error {
	psize := int64(C.NumElements() * A.Cols())
	Ar := A.FloatArray()
	ldA := A.LeadingIndex()
	Br := B.FloatArray()
	ldB := B.LeadingIndex()
	Cr := C.FloatArray()
	ldC := C.LeadingIndex()

	if nWorker <= 1 || psize <= limitOne {
		calgo.DMult(Cr, Ar, Br, alpha, beta, calgo.NOTRANS, ldC, ldA, ldB, B.Rows(),
			0, C.Cols(), 0, C.Rows(), vpLen, nB, mB)
		return nil
	}
	// here we have more than one worker available
	worker := func(cstart, cend, rstart, rend int, ready chan int) {
		calgo.DMult(Cr, Ar, Br, alpha, beta, calgo.NOTRANS, ldC, ldA, ldB, B.Rows(),
			cstart, cend, rstart, rend, vpLen, nB, mB)
		ready <- 1
	}
	colworks, rowworks := divideWork(C.Rows(), C.Cols(), nWorker)
	scheduleWork(colworks, rowworks, C.Cols(), C.Rows(), worker)
	return nil
}

// Calculate C = alpha*A.T*B + beta*C, C is M*N, A is P*M and B is P*N
func MMMultTransA(C, A, B *matrix.FloatMatrix, alpha, beta float64) error {
	psize := int64(C.NumElements() * B.Rows())
	Ar := A.FloatArray()
	ldA := A.LeadingIndex()
	Br := B.FloatArray()
	ldB := B.LeadingIndex()
	Cr := C.FloatArray()
	ldC := C.LeadingIndex()
	if nWorker <= 1 || psize <= limitOne {
		calgo.DMult(Cr, Ar, Br, alpha, beta, calgo.TRANSA, ldC, ldA, ldB,
			B.Rows(), 0, C.Cols(), 0, C.Rows(), vpLen, nB, mB)
		return nil
	}

	// here we have more than one worker available
	worker := func(cstart, cend, rstart, rend int, ready chan int) {
		calgo.DMult(Cr, Ar, Br, alpha, beta, calgo.TRANSA, ldC, ldA, ldB, B.Rows(),
			cstart, cend, rstart, rend, vpLen, nB, mB)
		ready <- 1
	}
	colworks, rowworks := divideWork(C.Rows(), C.Cols(), nWorker)
	scheduleWork(colworks, rowworks, C.Cols(), C.Rows(), worker)
	//scheduleWork(colworks, rowworks, worker)
	return nil
}

// Calculate C = alpha*A*B.T + beta*C, C is M*N, A is M*P and B is N*P
func MMMultTransB(C, A, B *matrix.FloatMatrix, alpha, beta float64) error {
	psize := int64(C.NumElements() * A.Cols())
	Ar := A.FloatArray()
	ldA := A.LeadingIndex()
	Br := B.FloatArray()
	ldB := B.LeadingIndex()
	Cr := C.FloatArray()
	ldC := C.LeadingIndex()
	if nWorker <= 1 || psize <= limitOne {
		calgo.DMult(Cr, Ar, Br, alpha, beta, calgo.TRANSB, ldC, ldA, ldB,
			B.Rows(), 0, C.Cols(), 0, C.Rows(), vpLen, nB, mB)
		return nil
	}

	// here we have more than one worker available
	worker := func(cstart, cend, rstart, rend int, ready chan int) {
		calgo.DMult(Cr, Ar, Br, alpha, beta, calgo.TRANSB, ldC, ldA, ldB, B.Rows(),
			cstart, cend, rstart, rend, vpLen, nB, mB)
		ready <- 1
	}
	colworks, rowworks := divideWork(C.Rows(), C.Cols(), nWorker)
	scheduleWork(colworks, rowworks, C.Cols(), C.Rows(), worker)
	//scheduleWork(colworks, rowworks, worker)
	return nil
}

// Calculate C = alpha*A.T*B.T + beta*C, C is M*N, A is P*M and B is N*P
func MMMultTransAB(C, A, B *matrix.FloatMatrix, alpha, beta float64) error {
	psize := int64(C.NumElements() * A.Rows())
	Ar := A.FloatArray()
	ldA := A.LeadingIndex()
	Br := B.FloatArray()
	ldB := B.LeadingIndex()
	Cr := C.FloatArray()
	ldC := C.LeadingIndex()
	if nWorker <= 1 || psize <= limitOne {
		calgo.DMult(Cr, Ar, Br, alpha, beta, calgo.TRANSA|calgo.TRANSB, ldC, ldA, ldB,
			B.Rows(), 0, C.Cols(), 0, C.Rows(), vpLen, nB, mB)
		return nil
	}

	// here we have more than one worker available
	worker := func(cstart, cend, rstart, rend int, ready chan int) {
		calgo.DMult(Cr, Ar, Br, alpha, beta, calgo.TRANSA|calgo.TRANSB, ldC, ldA, ldB,
			B.Rows(), cstart, cend, rstart, rend, vpLen, nB, mB)
		ready <- 1
	}
	colworks, rowworks := divideWork(C.Rows(), C.Cols(), nWorker)
	scheduleWork(colworks, rowworks, C.Cols(), C.Rows(), worker)
	//scheduleWork(colworks, rowworks, worker)
	return nil
}

// Calculate C = alpha*A*B + beta*C, C is M*N, A is M*M and B is M*N
func MMSymm(C, A, B *matrix.FloatMatrix, alpha, beta float64) error {

	if A.Rows() != A.Cols() {
		return errors.New("A matrix not square matrix.")
	}
	psize := int64(C.NumElements() * A.Cols())
	Ar := A.FloatArray()
	ldA := A.LeadingIndex()
	Br := B.FloatArray()
	ldB := B.LeadingIndex()
	Cr := C.FloatArray()
	ldC := C.LeadingIndex()

	if nWorker <= 1 || psize <= limitOne {
		calgo.DMultSymm(Cr, Ar, Br, alpha, beta, calgo.LEFT|calgo.LOWER, ldC, ldA, ldB,
			A.Cols(), 0, C.Cols(), 0, C.Rows(), vpLen, nB, mB)
		return nil
	}
	// here we have more than one worker available
	worker := func(cstart, cend, rstart, rend int, ready chan int) {
		calgo.DMultSymm(Cr, Ar, Br, alpha, beta, calgo.LEFT|calgo.LOWER, ldC, ldA, ldB,
			A.Cols(), cstart, cend, rstart, rend, vpLen, nB, mB)
		ready <- 1
	}
	colworks, rowworks := divideWork(C.Rows(), C.Cols(), nWorker)
	scheduleWork(colworks, rowworks, C.Cols(), C.Rows(), worker)
	return nil
}

// Calculate C = alpha*A*B + beta*C, C is M*N, A is M*M and B is M*N
func MMSymmUpper(C, A, B *matrix.FloatMatrix, alpha, beta float64) error {

	if A.Rows() != A.Cols() {
		return errors.New("A matrix not square matrix.")
	}
	psize := int64(C.NumElements() * A.Cols())
	Ar := A.FloatArray()
	ldA := A.LeadingIndex()
	Br := B.FloatArray()
	ldB := B.LeadingIndex()
	Cr := C.FloatArray()
	ldC := C.LeadingIndex()

	if nWorker <= 1 || psize <= limitOne {
		calgo.DMultSymm(Cr, Ar, Br, alpha, beta, calgo.LEFT|calgo.UPPER, ldC, ldA, ldB,
			A.Cols(), 0, C.Cols(), 0, C.Rows(), vpLen, nB, mB)
		return nil
	}
	// here we have more than one worker available
	worker := func(cstart, cend, rstart, rend int, ready chan int) {
		calgo.DMultSymm(Cr, Ar, Br, alpha, beta, calgo.LEFT|calgo.UPPER, ldC, ldA, ldB,
			A.Cols(), cstart, cend, rstart, rend, vpLen, nB, mB)
		ready <- 1
	}
	colworks, rowworks := divideWork(C.Rows(), C.Cols(), nWorker)
	scheduleWork(colworks, rowworks, C.Cols(), C.Rows(), worker)
	return nil
}

// Y = alpha*A.T*X + beta*Y
func MVMultTransA(Y, A, X *matrix.FloatMatrix, alpha, beta float64) error {

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
	lenY := Y.Rows()
	if Y.Rows() == 1 {
		// row vector
		incY = Y.LeadingIndex()
		lenY = Y.Cols()
	}
	Xr := X.FloatArray()
	incX := 1
	lenX := X.Rows()
	if X.Rows() == 1 {
		// row vector
		incX = X.LeadingIndex()
		lenX = X.Cols()
	}
	calgo.DMultMV(Yr, Ar, Xr, alpha, beta, calgo.TRANSA, incY, ldA, incX,
		0, lenX, 0, lenY, vpLen, mB)
	return nil
}

// A = A + alpha*X*Y.T; A is N*N symmetric, X is row or column vector of length N.
func MVSymmUpdateUpper(A, X *matrix.FloatMatrix, alpha float64) error {

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
	calgo.DSymmRankMV(Ar, Xr, alpha, calgo.UPPER, ldA, incX, 0, A.Cols(), 0)
	return nil
}

func MVSymm2UpdateUpper(A, X, Y *matrix.FloatMatrix, alpha float64) error {

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
	calgo.DSymmRank2MV(Ar, Xr, Yr, alpha, calgo.UPPER, ldA, incY, incX, 0, A.Cols(), 0)
	return nil
}
