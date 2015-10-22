// Copyright (c) Harri Rautila, 2012

// This file is part of go.opt/matrix package. It is free software, distributed
// under the terms of GNU Lesser General Public License Version 3, or any later
// version. See the COPYING tile included in this archive.

package matrix

import "math"

// Test for equality. Return true if for all i,j: all A[i,j] = B[i,j]
func (A *FloatMatrix) Equal(B *FloatMatrix) bool {
    if !A.SizeMatch(B.Size()) {
        return false
    }
    for k := 0; k < A.NumElements(); k++ {
        rka := realIndex(k, A.Rows(), A.LeadingIndex())
        rkb := realIndex(k, B.Rows(), B.LeadingIndex())
        if A.elements[rka] != B.elements[rkb] {
            return false
        }
    }
    return true
}

// Test for element wise less-than. Return true if for all i,j: A[i,j] < B[i,j]
func (A *FloatMatrix) Less(B *FloatMatrix) bool {
    if !A.SizeMatch(B.Size()) {
        return false
    }
    for k := 0; k < A.NumElements(); k++ {
        rka := realIndex(k, A.Rows(), A.LeadingIndex())
        rkb := realIndex(k, B.Rows(), B.LeadingIndex())
        if A.elements[rka] >= B.elements[rkb] {
            return false
        }
    }
    return true

}

// Test for element wise less-or-equal.
// Return true if for all i,j: A[i,j] <= B[i,j]
func (A *FloatMatrix) LessOrEqual(B *FloatMatrix) bool {
    if !A.SizeMatch(B.Size()) {
        return false
    }
    for k := 0; k < A.NumElements(); k++ {
        rka := realIndex(k, A.Rows(), A.LeadingIndex())
        rkb := realIndex(k, B.Rows(), B.LeadingIndex())
        if A.elements[rka] > B.elements[rkb] {
            return false
        }
    }
    return true
}

// Test for element wise greater-than.
// Return true if for all i,j: A[i,j] > B[i,j]
func (A *FloatMatrix) Greater(B *FloatMatrix) bool {
    if !A.SizeMatch(B.Size()) {
        return false
    }
    for k := 0; k < A.NumElements(); k++ {
        rka := realIndex(k, A.Rows(), A.LeadingIndex())
        rkb := realIndex(k, B.Rows(), B.LeadingIndex())
        if A.elements[rka] <= B.elements[rkb] {
            return false
        }
    }
    return true
}

// Test for element wise greater-than-or-equal.
// Return true if for all i,j: A[i,j] >= B[i,j]
func (A *FloatMatrix) GreaterOrEqual(B *FloatMatrix) bool {
    if !A.SizeMatch(B.Size()) {
        return false
    }
    for k := 0; k < A.NumElements(); k++ {
        rka := realIndex(k, A.Rows(), A.LeadingIndex())
        rkb := realIndex(k, B.Rows(), B.LeadingIndex())
        if A.elements[rka] < B.elements[rkb] {
            return false
        }
    }
    return true
}

// Default relative tolenrance, these are spyed from Python NumPy
const RTOL = 1.0000000000000001e-05

// Default absolute tolerance
const ATOL = 1e-8

// Return true if A is element-wise equal to with a tolerance. The tolerance
// values are positive, typically very small numbers. If tolerances parameter
// is not given default values are used. If one tolerance value is given, it is
// used as absolute tolerance. If two values are given, first is used as absolute
// tolerance and second as relative tolerance.
func (A *FloatMatrix) AllClose(B *FloatMatrix, tolerances ...float64) bool {
    if !A.SizeMatch(B.Size()) {
        return false
    }
    rtol := RTOL
    atol := ATOL
    if len(tolerances) == 1 {
        atol = tolerances[0]
    } else if len(tolerances) > 1 {
        atol = tolerances[0]
        rtol = tolerances[1]
    }
    for k := 0; k < A.NumElements(); k++ {
        rka := realIndex(k, A.Rows(), A.LeadingIndex())
        rkb := realIndex(k, B.Rows(), B.LeadingIndex())
        df := math.Abs(A.elements[rka] - B.elements[rkb])
        if df > atol+rtol*math.Abs(B.elements[rkb]) {
            return false
        }
    }
    return true
}

// Local Variables:
// tab-width: 4
// End:
