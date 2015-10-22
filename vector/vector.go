// Author: slowpoke <proxypoke at lavabit dot com>
// Repository: https://gist.github.com/proxypoke/vector
//
// This program is free software under the non-terms
// of the Anti-License. Do whatever the fuck you want.

// Package vector implements mathematical vectors over float64 values, with
// common operations defined on them, like addition, the scalar product or
// normalization. 
//
// Operations that result in a new vector (like addition or the
// cross product) have both an in-place and a non-destructive version, with the
// first being a method on the Vector type and the latter being a function.
//
// In-place operations on vectors return the vector to make it possible to
// chain ("pipe") operations. This is purely convenience.
package vector

import (
	"math"
)

// A vector over float64 values
type Vector struct {
	dims []float64 // the elements of the vector
	ndim uint      // the dimension of the vector
}

// ============================= [ Constructors ] =============================

// Create a Vector with dimension n, with all values initialized to 0.
func New(n uint) (v *Vector) {
	v = new(Vector)
	v.dims = make([]float64, n)
	v.ndim = n
	return
}

// Create a Vector from a slice. Its dimension is equal to len(slice).
func NewFrom(dims []float64) (v *Vector) {
	v = new(Vector)
	v.ndim = uint(len(dims))
	v.dims = make([]float64, v.ndim)
	copy(v.dims, dims)
	return
}

// Make a deep copy of the Vector.
func (v *Vector) Copy() *Vector {
	return NewFrom(v.dims)
}

// =========================== [ General Methods ] ============================

// Get the dimension of the Vector.
func (v Vector) Dim() uint {
	return v.ndim
}

// Get the value of the nth element in the Vector.
func (v Vector) Get(n uint) (val float64, err error) {
	if n >= v.Dim() {
		err = IndexError(n)
		return
	}
	val = v.dims[n]
	return
}

// Set the value of the nth element in the Vector.
func (v Vector) Set(n uint, x float64) (err error) {
	if n >= v.Dim() {
		err = IndexError(n)
		return
	}
	v.dims[n] = x
	return
}

// Calculate the length of the Vector.
func (v Vector) Len() (result float64) {
	for _, val := range v.dims {
		result += math.Pow(val, 2)
	}
	result = math.Sqrt(result)
	// Account for floating point imprecison
	// XXX: This is probably a bad solution, but it works for now.
	ε := 1.00000000000000000005
	if math.Abs(1-result) < ε {
		result = 1
	}
	return
}

// ========================= [ In-place operations ] ==========================

// Add another Vector, in-place.
func (v *Vector) Add(other *Vector) (ret *Vector, err error) {
	err = checkDims(v, other)
	if err == nil {
		for i := range v.dims {
			v.dims[i] += other.dims[i]
		}
		ret = v
	}
	return
}

// Substract another Vector, in-place.
func (v *Vector) Substract(other *Vector) (ret *Vector, err error) {
	err = checkDims(v, other)
	if err == nil {
		for i := range v.dims {
			v.dims[i] -= other.dims[i]
		}
		ret = v
	}
	return
}

// In-place scalar multiplication.
func (v *Vector) Scale(x float64) *Vector {
	for i := range v.dims {
		v.dims[i] *= x
	}
	return v
}

// Normalize the Vector (length == 1). In-place.
func (v *Vector) Normalize() *Vector {
	l := v.Len()
	for i := range v.dims {
		v.dims[i] /= l
	}
	return v
}

// Cross product with another vector, in-place.
// Returns error when ndim of either vector != 3.
func (v *Vector) CrossProduct(other *Vector) (*Vector, error) {
	if v.ndim != 3 || other.ndim != 3 {
		err := CrossError{v.ndim, other.ndim}
		return nil, err
	}
	x := v.dims[1]*other.dims[2] - v.dims[2]*other.dims[1]
	y := v.dims[2]*other.dims[0] - v.dims[0]*other.dims[2]
	z := v.dims[0]*other.dims[1] - v.dims[1]*other.dims[0]
	v.dims[0] = x
	v.dims[1] = y
	v.dims[2] = z
	return v, nil
}

// ============================== [ Functions ] ===============================

// Compare two vectors. Returns true if all values are the same.
func Equal(a, b *Vector) (equal bool) {
	err := checkDims(a, b)
	if err == nil {
		equal = true
		for i := range a.dims {
			if a.dims[i] != b.dims[i] {
				equal = false
			}
		}
	}
	return
}

// Add two Vectors, returning a new Vector.
func Add(a, b *Vector) (*Vector, error) {
	return a.Copy().Add(b)
}

// Substract two Vectors, returning new Vector.
func Substract(a, b *Vector) (*Vector, error) {
	return a.Copy().Substract(b)
}

// Scalar multiplication of a Vector, returning a new Vector.
func Scale(v *Vector, x float64) *Vector {
	return v.Copy().Scale(x)
}

// Normalize a vector, returning a new Vector.
func Normalize(v *Vector) *Vector {
	return v.Copy().Normalize()
}

// Dot-product of two Vectors.
func DotProduct(a, b *Vector) (dot float64, err error) {
	for i := range a.dims {
		dot += a.dims[i] * b.dims[i]
	}
	return
}

// Angle Θ (theta) between two vectors.
func Angle(a, b *Vector) (Θ float64, err error) {
	err = checkDims(a, b)
	if err == nil {
		norm_a := Normalize(a)
		norm_b := Normalize(b)
		dot, _ := DotProduct(norm_a, norm_b)
		Θ = math.Acos(dot)
	}
	return
}

// Cross product of two vectors.
// Returns error when ndim of either vector != 3.
func CrossProduct(a, b *Vector) (*Vector, error) {
	return a.Copy().CrossProduct(b)
}

// =========================== [ Helper Functions ] ===========================

// Check if two vectors have the same dimension.
func checkDims(a, b *Vector) (err error) {
	if a.ndim != b.ndim {
		err = DimError{a.ndim, b.ndim}
	}
	return
}
