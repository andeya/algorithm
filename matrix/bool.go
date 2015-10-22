// Copyright (c) Harri Rautila, 2012

// This file is part of go.opt/matrix package. It is free software, distributed
// under the terms of GNU Lesser General Public License Version 3, or any later
// version. See the COPYING tile included in this archive.

package matrix

// Test if matrices A and B are equal. Returns false if sizes of the matrices
// do not match or if any A[i,j] != B[i,j]
func Equal(A, B Matrix) bool {
    if !A.SizeMatch(B.Size()) {
        return false
    }
    if !EqualTypes(A, B) {
        return false
    }
    switch A.(type) {
    case *FloatMatrix:
        Ar := A.(*FloatMatrix).FloatArray()
        Br := B.(*FloatMatrix).FloatArray()
        for i, v := range Ar {
            if Br[i] != v {
                return false
            }
        }
    case *ComplexMatrix:
        Ar := A.(*ComplexMatrix).FloatArray()
        Br := B.(*ComplexMatrix).FloatArray()
        for i, v := range Ar {
            if Br[i] != v {
                return false
            }
        }
    }
    return true
}

// Test if matrices A and B are not equal. Returns true if matrices are of different size
// or if not all element values are equal.
func NotEqual(A, B Matrix) bool {
    if !A.SizeMatch(B.Size()) {
        return true
    }
    return !Equal(A, B)
}

// Test if matrices are of same type. Return true if all are same type
// as the first element, otherwise return false. For parameter count <=1
// return true.
func EqualTypes(As ...Matrix) bool {
    if len(As) <= 1 {
        return true
    }
    return As[0].EqualTypes(As...)
}

// Local Variables:
// tab-width: 4
// End:
