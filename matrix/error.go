// Copyright (c) Harri Rautila, 2012

// This file is part of go.opt/matrix package. It is free software, distributed
// under the terms of GNU Lesser General Public License Version 3, or any later
// version. See the COPYING tile included in this archive.

package matrix

import (
    "errors"
    "fmt"
)

const (
    //The matrix returned was nil.
    errorNilMatrix = iota + 1
    //The dimensions of the inputs do not make sense for this operation.
    errorDimensionMismatch
    //The indices provided are out of bounds.
    errorIllegalIndex
    //The matrix provided has a singularity.
    exceptionSingular
    //The matrix provided is not positive semi-definite.
    exceptionNotSPD
    // The argument type mismatch
    exceptionIllegalType
)

type error_ int

func (e error_) Error() string {
    switch e {
    case errorNilMatrix:
        return "Matrix is nil"
    case errorDimensionMismatch:
        return "Input dimensions do not match"
    case errorIllegalIndex:
        return "Index out of bounds"
    case exceptionSingular:
        return "Matrix is singular"
    case exceptionNotSPD:
        return "Matrix is not positive semidefinite"
    case exceptionIllegalType:
        return "Matrix type is not appropriate"
    }
    return fmt.Sprintf("Unknown error code %d", e)
}
func (e error_) String() string {
    return e.Error()
}

var (
    //The matrix returned was nil.
    ErrorNilMatrix error_ = error_(errorNilMatrix)
    //The dimensions of the inputs do not make sense for this operation.
    ErrorDimensionMismatch error_ = error_(errorDimensionMismatch)
    //The indices provided are out of bounds.
    ErrorIllegalIndex error_ = error_(errorIllegalIndex)
    //The matrix provided has a singularity.
    ExceptionSingular error_ = error_(exceptionSingular)
    //The matrix provided is not positive semi-definite.
    ExceptionNotSPD error_ = error_(exceptionNotSPD)
    //Not corrent type
    ExceptionIllegalType error_ = error_(exceptionIllegalType)
)

var panicOnError bool = false

func PanicOnError(flag bool) {
    panicOnError = flag
}

func onError(msg string) error {
    if panicOnError {
        panic(msg)
    }
    return errors.New(msg)
}

// Local Variables:
// tab-width: 4
// indent-tabs-mode: nil
// End:
