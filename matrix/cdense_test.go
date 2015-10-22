// Copyright (c) Harri Rautila, 2012

// This file is part of go.opt/matrix package. It is free software, distributed
// under the terms of GNU Lesser General Public License Version 3, or any later
// version. See the COPYING tile included in this archive.

package matrix

import (
    //"fmt"
    "testing"
)

func TestA(t *testing.T) {
    t.Logf("Test complex matrix printing.\n")
    A := ComplexZeros(2, 2)
    t.Logf("A:\n%v\n", A)
}

func TestCParse(t *testing.T) {
    t.Logf("Test matrix string parsing.\n")
    s := `[(1.0+0i) (+2-1i) (3.0+0i); ( 4.2-.5i) (-5 -.1i) (6+0i)]`
    A, err := ComplexParse(s)
    if err != nil {
        t.Fatal(err)
    }
    t.Logf("A :\n%v\n", A)
    t.Logf("A size: %d rows, %d cols\n", A.Rows(), A.Cols())
    t.Logf("A elems: %d\n", A.NumElements())
    D := A.Transpose()
    t.Logf("D = A.transpose:\n%v\n", D)
    r, c := D.Size()
    t.Logf("D size: %d rows, %d cols\n", r, c)
}

func TestRand(t *testing.T) {
    t.Logf("Test matrix creation.\n")
    A := ComplexUniform(3, 2)
    t.Logf("A :\n%v\n", A)
    B := ComplexUniformSymmetric(2)
    t.Logf("B :\n%v\n", B)
}

/*
func TestCopy(t *testing.T) {
	t.Logf("Test creating and setting elements.\n")
	A := MakeComplexMatrix(2, 3, []float64{1,4,2,5,3,6})
	t.Logf("A:\n%v\n", A)
	C := A.Copy()
	t.Logf("C:\n%v\n", C)
	C.Set(0, 1, C.Get(0, 1)*10)
	B := MakeComplexMatrix(3, 2, []float64{1,2,3,4,5,6})
	t.Logf("B:\n%v\n", B)
}

func TestBool(t *testing.T) {
	t.Logf("Test matrix boolean operators.\n")
	A := MakeComplexMatrix(3, 2, []float64{1,4,2,5,3,6})
	t.Logf("A:\n%v\n", A)
	B := A.Copy()
	//B := MakeComplexMatrix(3, 2, []float64{1,2,3,4,5,6})
	t.Logf("B:\n%v\n", B)
	t.Logf("A == B: %v\n", A.Equal(B))
	t.Logf("A != B: %v\n", A.NotEqual(B))
}

func TestMath(t *testing.T) {
	t.Logf("Test matrix basic math.\n")
	A := Zeros(2, 2)
	t.Logf("A\n%v\n", A)
	A.Add(1.0)
	t.Logf("A += 1.0\n%v\n", A)
	A.Mult(9.0)
	t.Logf("A *= 9.0\n%v\n", A)
	A.Sub(1.0)
	t.Logf("A -= 1.0\n%v\n", A)
	A.Div(2.0)
	t.Logf("A /= 2.0\n%v\n", A)
	A.Remainder(3.0)
	t.Logf("A %%= 3.0\n%v\n", A)
	A.Neg()
	t.Logf("A = -A:\n%v\n", A)
	C := A.Times(A)
	t.Logf("C = A*A:\n%v\n", C)
	D := C.Plus(A)
	t.Logf("D = C+A:\n%v\n", D)
	F := D.Minus(A)
	t.Logf("F = D-A:\n%v\n", F)
	G := Zeros(3, 2); G.Add(1.0)
	H := G.Copy().Transpose()
	t.Logf("G:\n%v\n", G)
	t.Logf("H:\n%v\n", H)
	K := G.Times(H)
	t.Logf("K = G*H:\n%v\n", K)
}

func TestFuncs(t *testing.T) {
	t.Logf("Test matrix element wise math.\n")
	A := Zeros(2, 3)
	AddTwo := func (n float64) float64 {
		return n+2.0
	}
	C := A.Apply(AddTwo)
	t.Logf("C = AddTwo(A):\n%v\n", C)
}

func TestFuncs2(t *testing.T) {
	A := Ones(2, 3)
	B := Exp(A)
	t.Logf("B = Exp(A):\n%v\n", B)
}
*/

// Local Variables:
// tab-width: 4
// End:
