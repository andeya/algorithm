package lm

import "testing"
import "fmt"
import "math"
import "math/rand"

func TestWls(t *testing.T) {
	const n, p = 10, 2
	methods := [...]uint8{'q', 'c'}
	// Testing with y = x
	X := make([]float64, n*p)
	y := make([]float64, n)
	w := make([]float64, n)

	// Fill X matrix and simulate y's
	for i := 0; i < n; i++ {
		X[i] = 1.
		X[i+n] = (float64)(i)
		y[i] = X[i+n]
		w[i] = 1.
	}

	for _, method := range methods {
		// Run regression
		coef, status := Wls(X, n, p, y, w, method)

		// Print results
		fmt.Printf("Status: %d\n", status)
		fmt.Printf("Beta hat: %v\n", coef)

		// Check numerical accuracy
		l2Error := math.Sqrt(math.Pow(coef[0], 2) + math.Pow(coef[1]-1, 2))
		maxError := math.Sqrt(2) * math.Sqrt(math.Nextafter(1., 2.)-1.)
		if status > 0 || l2Error > maxError {
			t.Errorf("Method %c\tStatus %d\tL2 error %g > %g",
				method, status, l2Error, maxError)
		} else {
			fmt.Printf("Method %c\tL2 error %g < %g\n",
				method, l2Error, maxError)
		}
	}
}

func TestLmT(t *testing.T) {
	const n, p = 100, 2
	// Testing with y = x + e; x = 0, ..., n; e ~ t_2
	X := make([]float64, n*p)
	y := make([]float64, n)

	// Fill X matrix and simulate y's
	for i := 0; i < n; i++ {
		X[i] = 1.
		X[i+n] = (float64)(i)
		y[i] = X[i+n] + rand.NormFloat64()/math.Sqrt(rand.ExpFloat64())
	}

	// Run regression
	coef, tau, iterations, logLikelihood := LmT(X, n, p, y, 2., 100, 1e-8, 'q')

	// Print results
	fmt.Printf("Beta hat: %v\n", coef)
	fmt.Printf("Tau hat: %v\n", tau)
	fmt.Printf("Iterations: %v\n", iterations)
	fmt.Printf("Log likelihood: %v\n", logLikelihood)
}
