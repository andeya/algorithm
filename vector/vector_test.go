// Author: slowpoke <proxypoke at lavabit dot com>
// Repository: https://gist.github.com/proxypoke/vector
//
// This program is free software under the non-terms
// of the Anti-License. Do whatever the fuck you want.

package vector

import (
	"math"
	"math/rand"
	"testing"
)

// ========================== [ Constructor Tests ] ===========================

// Creates vectors with dimension from 0 to 99, checks if they actually have
// that dimension, then checks if the values are correctly initialized to 0.
func TestNew(t *testing.T) {
	var i, j uint
	for i = 0; i < 100; i++ {
		v := New(i)
		if v.Dim() != i {
			t.Errorf("Wrong dimension. Got %d, expected %d.", v.Dim(), i)
		}
		for j = 0; j < i; j++ {
			// XXX: If the Get method errors, this test will still pass. This
			// is because Get() would then return an uninitialized float64 for
			// val, which is 0 and therefore what the test expects.
			val, _ := v.Get(j)
			if val != 0 {
				t.Error("Newly initialized vector has a value != 0.")
			}
		}
	}
}

// Creates vectors with randomized slices, then checks whether they have the
// correct dimension (len(slice)) and whether they have been correctly
// initialized.
func TestNewFrom(t *testing.T) {
	var i, j uint
	for i = 0; i < 100; i++ {
		randslice := makeRandSlice(i)
		v := NewFrom(randslice)
		if v.Dim() != i {
			t.Errorf("Wrong dimension. Got %d, expected %d.", v.Dim(), i)
		}
		for j = 0; j < i; j++ {
			val, _ := v.Get(j)
			if val != randslice[j] {
				t.Error(
					"Wrong values in vector initialized from a random slice.")
			}
		}
	}
}

// Creates pseudo-random vectors with various dimensions, copies them and
// verifies that the new vector is equal.
func TestCopy(t *testing.T) {
	var i uint
	for i = 0; i < 100; i++ {
		v := makeRandomVector(i)
		w := v.Copy()
		if !Equal(v, w) {
			t.Error("Copied vector is not equal to source vector.")
		}
	}
}

// =================== [ General Methods/Functions Tests ] ====================

// Creates pseudo-random vectors with various dimensions, then check if Get()
// returns the correct values and errors on out-of-range indexes.
func TestGet(t *testing.T) {
	var i uint
	for i = 0; i < 100; i++ {
		v := makeRandomVector(i)
		for j, val := range v.dims {
			getval, err := v.Get(uint(j))
			if err != nil {
				t.Error("Get() errored on a correct index.")
			}
			if val != getval {
				t.Error("Get() returned a wrong value.")
			}
		}
		_, err := v.Get(v.Dim())
		if err == nil {
			t.Error("Get didn't error on an out-of-range index.")
		}
	}
}

// Creates uninitialized vectors of various dimensions, then sets their values
// to pseudo-random values. It then compares those values to check if they
// were set correctly. Also verifies is Set() correctly errors on out-of-range
// indexes.
func TestSet(t *testing.T) {
	var i, j uint
	for i = 0; i < 100; i++ {
		v := New(i)
		for j = 0; j < i; j++ {
			val := rand.ExpFloat64()
			err := v.Set(j, val)
			if err != nil {
				t.Error("Set() errored on a correct index.")
			}
			if v.dims[j] != val {
				t.Error("Set didn't correctly set a value.")
			}
		}
		err := v.Set(v.Dim(), 0)
		if err == nil {
			t.Error("Set didn't error on an out-of-range index.")
		}
	}
}

// Creates a vector with known length, then compares the expected value with
// what Len() returns.
func TestLen(t *testing.T) {
	v := New(1)
	v.Set(0, 2) // has length 2
	if v.Len() != 2 {
		t.Error("Len returned a wrong length")
	}
}

// Creates Vectors which are known to be (un)equal, then verifies that Equal()
// has correct oytput.
func TestEqual(t *testing.T) {
	slc := make([]float64, 10)
	for i := range slc {
		slc[i] = float64(i)
	}

	v := NewFrom(slc)
	w := NewFrom(slc)
	if !Equal(v, w) {
		t.Error("Equal() != true for equal vectors.")
	}

	w = New(10)
	if Equal(v, w) {
		t.Error("Equal() == true for unequal vectors.")
	}
}

// =========================== [ Operation Tests ] ============================

// Creates pesudo-random vectors, then adds them first as a non-destructive,
// then as an in-place operations, checking if both operation were correct.
func TestAdd(t *testing.T) {
	var i, j uint
	for i = 1; i < 100; i++ {
		a := makeRandomVector(i)
		b := makeRandomVector(i)
		c, _ := Add(a, b)

		for j = 0; j < i; j++ {
			if c.dims[j] != a.dims[j]+b.dims[j] {
				t.Error("Addition failed, didn't get expected values.")
				t.Logf("%f + %f != %f", a.dims[j], b.dims[j], c.dims[j])
			}
		}

		// Test in-place addition.
		c = a.Copy()
		c.Add(b)

		for j = 0; j < i; j++ {
			if c.dims[j] != a.dims[j]+b.dims[j] {
				t.Error(
					"In-place Addition failed, didn't get expected values.")
				t.Logf("%f + %f != %f", a.dims[j], b.dims[j], c.dims[j])
			}
		}
	}
}

// Same as TestAdd, but with substraction. Heck, it's basically the same code.
func TestSubstract(t *testing.T) {
	var i, j uint
	for i = 1; i < 100; i++ {
		a := makeRandomVector(i)
		b := makeRandomVector(i)
		c, _ := Substract(a, b)

		for j = 0; j < i; j++ {
			if c.dims[j] != a.dims[j]-b.dims[j] {
				t.Error("Substraction failed, didn't get expected values.")
				t.Logf("%f - %f != %f", a.dims[j], b.dims[j], c.dims[j])
			}
		}

		// Test in-place sybstraction
		c = a.Copy()
		c.Substract(b)

		for j = 0; j < i; j++ {
			if c.dims[j] != a.dims[j]-b.dims[j] {
				t.Error(
					"In-place Substraction failed, didn't get expected values.")
				t.Logf("%f - %f != %f", a.dims[j], b.dims[j], c.dims[j])
			}
		}
	}
}

// Creates pseudo-random vectors, does scalar multiplication with pseudo-random
// floats, then checks if the result is correct. It checks both the in-place
// and the non-destructive version.
func TestScale(t *testing.T) {
	var i, j uint
	for i = 0; i < 100; i++ {
		a := makeRandomVector(i)
		x := rand.ExpFloat64()
		b := Scale(a, x)

		for j = 0; j < i; j++ {
			if b.dims[j] != a.dims[j]*x {
				t.Error("Scalar Multiplication failed, ",
					"didn't get expected values.")
				t.Logf("%f * %f != %f", a.dims[j], x, b.dims[j])
			}
		}

		// Test in-place scalar multiplication
		b = a.Copy()
		b.Scale(x)

		for j = 0; j < i; j++ {
			if b.dims[j] != a.dims[j]*x {
				t.Error("In-place Scalar Multiplication failed, ",
					"didn't get expected values.")
				t.Logf("%f * %f != %f", a.dims[j], x, b.dims[j])
			}
		}
	}
}

// Creates pseudo-random vectors, normalizes them both in-place and
// non-destructive, and verifies that the result is correct.
func TestNormalize(t *testing.T) {
	var i uint
	// It makes no sense to normalize a zero vector, therefore we start at 1.
	for i = 1; i < 100; i++ {
		a := makeRandomVector(i)
		b := Normalize(a)

		if b.Len() != float64(1) {
			t.Error("Normalization failed, vector doesn't have length 1.")
			t.Logf("%f != 1", b.Len())
		}
	}
}

// Uses vectors with known angles to calculate their DotProduct, then verifies
// if the result is correct.
func TestDotProduct(t *testing.T) {
	a := New(2)
	b := New(2)

	// Set the vectors as parallel.
	a.Set(0, 1)
	b.Set(0, 1)
	dot, _ := DotProduct(a, b)
	if dot != 1 {
		t.Error("Dot Product of parallel vectors isn't 1.")
	}

	// Set the vectors as orthogonal.
	b = New(2)
	b.Set(1, 1)
	dot, _ = DotProduct(a, b)
	if dot != 0 {
		t.Error("Dot Product of orthogonal vectors isn't 0.")
	}

	// Set the vectors as anti-parallel.
	b = New(2)
	b.Set(0, -1)
	dot, _ = DotProduct(a, b)
	if dot != -1 {
		t.Error("Dot Product of anti-parallel vectors isn't -1.")
	}
}

// Uses vectors with known angles to verify that Angle() is correct.
func TestAngle(t *testing.T) {
	a := New(2)
	b := New(2)

	// Set the vectors as parallel (Θ == 0).
	a.Set(0, 1)
	b.Set(0, 1)
	Θ, _ := Angle(a, b)
	if Θ != 0 {
		t.Error("Angle between parallel vectors isn't 0.")
		t.Logf("%f != 0", Θ)
	}

	// Set the vectors as orthogonal (Θ == 0.5π).
	b = New(2)
	b.Set(1, 1)
	Θ, _ = Angle(a, b)
	if Θ != 0.5*math.Pi {
		t.Error("Angle between orthonal vectors isn't 0.5π.")
		t.Logf("%f != %f", Θ, 0.5*math.Pi)
	}

	// Set the vectors as anti-parallel (Θ == π).
	b = New(2)
	b.Set(0, -1)
	Θ, _ = Angle(a, b)
	if Θ != math.Pi {
		t.Error("Angle between anti-parallel vectors isn't π.")
		t.Logf("%f != %f", Θ, math.Pi)
	}
}

// Calculates the cross product of two pseudo-random vectors, then checks if
// the resulting vector is orthogonal to both the original vectors. Tests both
// in-place and non-destructive versions of the operation.
func TestCrossProduct(t *testing.T) {
	check := func(a, b, c *Vector) {
		dot_a, _ := DotProduct(a, c)
		dot_b, _ := DotProduct(b, c)
		ε := 0.0000000005
		if math.Abs(0-dot_a) < ε {
			dot_a = 0
		}
		if math.Abs(0-dot_b) < ε {
			dot_b = 0
		}
		if dot_a != 0 || dot_b != 0 {
			t.Error("Either or both vectors aren't orthogonal",
				"to their Cross Product.")
			t.Logf("a * c = %f", dot_a)
			t.Logf("b * c = %f", dot_b)
		}
	}

	a := makeRandomVector(3)
	b := makeRandomVector(3)
	c, _ := CrossProduct(a, b)

	check(a, b, c)

	// Check in-place, too.
	c = a.Copy()
	c.CrossProduct(b)

	check(a, b, c)

	// Check if vectors ∉ ℝ³ are rejected.
	d := New(2)
	e := New(4)
	_, err := CrossProduct(d, e)
	if err == nil {
		t.Error("CrossProduct() didn't error with invalid input vectors",
			"(∉ ℝ³)")
	}
}

// Check whether the various functions that take more than one vector error on
// being supplied with vectors of missmatched dimensions.
// It suffices to check the helper function checkDims, since every function
// must call it to verify its inputs.
func TestMissmatchedDims(t *testing.T) {
	a := New(2)
	b := New(3)

	err := checkDims(a, b)
	if err == nil {
		t.Error("Missmatched dimension check succeeded on unequal dimensions.")
	}

	a = New(4)
	b = New(4)
	err = checkDims(a, b)
	if err != nil {
		t.Error("Missmatched dimension check failed on equal dimensions.")
	}
}

// =========================== [ Helper Functions ] ===========================

// Helper function, makes pseudo-random slices.
func makeRandSlice(length uint) (randslice []float64) {
	randslice = make([]float64, length)
	for i := range randslice {
		randslice[i] = rand.ExpFloat64()
	}
	return
}

// Helper function, make a pseudo-random Vector with dimension dim.
func makeRandomVector(dim uint) *Vector {
	return NewFrom(makeRandSlice(dim))
}
