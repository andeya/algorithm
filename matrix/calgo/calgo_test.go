
package calgo

import (
    "testing"
)

const NROWS = 2000
var A, B, C []float64

func TestMakeData(t *testing.T) {
    // all values set matrix
    A = make([]float64, NROWS*NROWS)
    for i, _ := range A {
        A[i] = 2.0
    }
    // diagonal matrix
    B = make([]float64, NROWS*NROWS)
    for i := 0; i < NROWS; i++ {
        B[i*NROWS+i] = 1.0
    }
    C = make([]float64, NROWS*NROWS)
    t.Logf("A [%d,%d] non-zero matrix\n", NROWS, NROWS)
    t.Logf("B [%d,%d] diagonal matrix\n", NROWS, NROWS)
    t.Logf("C [%d,%d] result matrix\n", NROWS, NROWS)
    
}

func TestMultAB(t *testing.T) {
    // this should run quickly
    Mult(C, A, B, 1.0, NROWS, NROWS, NROWS, NROWS, NROWS, NROWS)
}

func TestMultBA(t *testing.T) {
    // this take more time
    Mult(C, B, A, 1.0, NROWS, NROWS, NROWS, NROWS, NROWS, NROWS)
}

// Local Variables:
// tab-width: 4
// indent-tabs-mode: nil
// End:
