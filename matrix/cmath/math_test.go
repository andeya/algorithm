package cmath

import (
	"fmt"
	"github.com/henrylee2cn/analysis/matrix"
	"testing"
)

func TestFMath(t *testing.T) {
	fmt.Printf("Test matrix basic math.\n")
	A := matrix.ComplexZeros(2, 2)
	fmt.Printf("A\n%v\n", A)
	A = Add(A, 1.0)
	fmt.Printf("A += 1.0\n%v\n", A)
	A = Scale(A, 9.0)
	fmt.Printf("A *= 9.0\n%v\n", A)
	A = Add(A, -1.0)
	fmt.Printf("A -= 1.0\n%v\n", A)
	A = Mul(A, A)
	fmt.Printf("A = A .* A\n%v\n", A)
}

// Local Variables:
// tab-width: 4
// End:
