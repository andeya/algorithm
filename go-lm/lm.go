package lm

// #cgo CFLAGS: -O2 -Wall
// #cgo LDFLAGS: -lf77blas -llapack -latlas -lm -lgfortran
// #include "wls.h"
// #include "lmT.h"
import "C"
import "unsafe"

// X must be in column-major order
func Wls(X []float64, n int, p int, y []float64, w []float64, method uint8) (
	[]float64, int) {
	// Allocate memory for intermediate objects
	XTX := make([]float64, p*p)
	sqw := make([]float64, n)
	sqwX := make([]float64, n*p)
	sqwy := make([]float64, n)
	coef := make([]float64, p)

	// Call C function for WLS fit
	status := C.wls(
		(*C.double)(unsafe.Pointer(&X[0])), C.int(n), C.int(p),
		(*C.double)(unsafe.Pointer(&y[0])),
		(*C.double)(unsafe.Pointer(&w[0])),
		C.char(method),
		(*C.double)(unsafe.Pointer(&XTX[0])),
		(*C.double)(unsafe.Pointer(&sqw[0])),
		(*C.double)(unsafe.Pointer(&sqwX[0])),
		(*C.double)(unsafe.Pointer(&sqwy[0])),
		(*C.double)(unsafe.Pointer(&coef[0])))

	return coef, (int)(status)
}

func LmT(X []float64, n int, p int, y []float64,
	nu float64, maxIter int, tol float64, method uint8) (
	coef []float64, tau float64, iterations int, logLikelihood float64) {
	// Setup variables for results
	coef = make([]float64, p)

	// Call C code
	iterations = (int)(C.lmT(
		(*C.double)(unsafe.Pointer(&X[0])),
		C.int(n), C.int(p),
		(*C.double)(unsafe.Pointer(&y[0])),
		C.double(nu),
		C.int(maxIter), C.double(tol), C.char(method),
		(*C.double)(unsafe.Pointer(&logLikelihood)),
		(*C.double)(unsafe.Pointer(&coef[0])),
		(*C.double)(unsafe.Pointer(&tau))))

	// Return results
	return
}
