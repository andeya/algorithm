// Copyright (c) Harri Rautila, 2012

// This file is part of github.com/hrautila/linalg package. 
// It is free software, distributed under the terms of GNU Lesser General Public 
// License Version 3, or any later version. See the COPYING tile included in this archive.

package linalg

import (
    "fmt"
    "strings"
)

// Type Opt holds one BLAS/LAPACK index or parameter option.
type Opt struct {
    Name string
    Val  int
}

// LinalgIndex structure holds fields for various BLAS/LAPACK indexing
// variables.
type IndexOpts struct {
    // these for BLAS and LAPACK
    N       int // default: -1
    Nx      int // default: -1
    Ny      int // default: -1
    M       int // default: -1
    Ma      int // default: -1
    Mb      int // default: -1
    LDa     int // default: 0
    LDb     int // default: 0
    LDc     int // default: 0
    IncX    int // default: 0
    IncY    int // default: 0
    OffsetX int // default: 0
    OffsetY int // default: 0
    OffsetA int // default: 0
    OffsetB int // default: 0
    OffsetC int // default: 0
    K       int // default: -1
    Ku      int // default: -1
    Kl      int // default: 0
    // these used in LAPACK
    Nrhs     int // default: -1
    OffsetD  int // default: 0
    OffsetDL int // default: 0
    OffsetDU int // default: 0
    LDw      int // default: 0
    LDz      int // default: 0
    OffsetW  int // default: 0
    OffsetZ  int // default: 0
    LDu      int // default: 0
    LDvt     int // default: 0
    LDt      int // default: 0
    OffsetS  int // default: 0
    OffsetU  int // default: 0
    OffsetVt int // default: 0
}

// Parse option list and return index structure with relevant fields set and
// other fields with default values.
func GetIndexOpts(opts ...Option) *IndexOpts {
    is := &IndexOpts{
        -1, -1, -1, // n, nX, nY
        -1, -1, -1, // m, mA, mB
        0, 0, 0, // ldA, ldB, ldC
        1, 1, // incX, incY
        0, 0, 0, 0, 0, // offsetX, ... offsetC
        -1, -1, 0, // k, ku, kl
        -1,      // nrhs
        0, 0, 0, // offsetD, offsetDL, OffsetDU,
        0, 0, // LDw, LDz
        0, 0, // OffsetW, OffsetZ
        0, 0, 0, // LDu, LDvt, LDt
        0, 0, 0, // OffsetS, OffsetU, OffsetVt
    }

loop:
    for _, o := range opts {
        if _, ok := o.(*IOpt); !ok {
            continue loop
        }
        switch {
        case strings.EqualFold(o.Name(), "inc"):
            is.IncX = o.Int()
            is.IncY = o.Int()
        case strings.EqualFold(o.Name(), "incx"):
            is.IncX = o.Int()
        case strings.EqualFold(o.Name(), "incy"):
            is.IncY = o.Int()
        case strings.EqualFold(o.Name(), "lda"):
            is.LDa = o.Int()
        case strings.EqualFold(o.Name(), "ldb"):
            is.LDb = o.Int()
        case strings.EqualFold(o.Name(), "ldc"):
            is.LDc = o.Int()
        case strings.EqualFold(o.Name(), "ldw"):
            is.LDw = o.Int()
        case strings.EqualFold(o.Name(), "ldz"):
            is.LDz = o.Int()
        case strings.EqualFold(o.Name(), "ldu"):
            is.LDu = o.Int()
        case strings.EqualFold(o.Name(), "ldvt"):
            is.LDvt = o.Int()
        case strings.EqualFold(o.Name(), "ldt"):
            is.LDt = o.Int()
        case strings.EqualFold(o.Name(), "offset"):
            is.OffsetX = o.Int()
            is.OffsetY = o.Int()
            is.OffsetA = o.Int()
            is.OffsetB = o.Int()
            is.OffsetC = o.Int()
        case strings.EqualFold(o.Name(), "offsetx"):
            is.OffsetX = o.Int()
        case strings.EqualFold(o.Name(), "offsety"):
            is.OffsetY = o.Int()
        case strings.EqualFold(o.Name(), "offseta"):
            is.OffsetA = o.Int()
        case strings.EqualFold(o.Name(), "offsetb"):
            is.OffsetB = o.Int()
        case strings.EqualFold(o.Name(), "offsetc"):
            is.OffsetC = o.Int()
        case strings.EqualFold(o.Name(), "offsetw"):
            is.OffsetW = o.Int()
        case strings.EqualFold(o.Name(), "offsetd"):
            is.OffsetD = o.Int()
        case strings.EqualFold(o.Name(), "offsetdl"):
            is.OffsetDL = o.Int()
        case strings.EqualFold(o.Name(), "offsetdu"):
            is.OffsetDU = o.Int()
        case strings.EqualFold(o.Name(), "offsetdw"):
            is.OffsetW = o.Int()
        case strings.EqualFold(o.Name(), "offsetdz"):
            is.OffsetZ = o.Int()
        case strings.EqualFold(o.Name(), "offsetu"):
            is.OffsetU = o.Int()
        case strings.EqualFold(o.Name(), "offsets"):
            is.OffsetS = o.Int()
        case strings.EqualFold(o.Name(), "offsetvt"):
            is.OffsetVt = o.Int()
        case strings.EqualFold(o.Name(), "n"):
            is.N = o.Int()
            is.Nx = o.Int()
            is.Ny = o.Int()
        case strings.EqualFold(o.Name(), "nx"):
            is.Nx = o.Int()
        case strings.EqualFold(o.Name(), "ny"):
            is.Ny = o.Int()
        case strings.EqualFold(o.Name(), "m"):
            is.M = o.Int()
            is.Ma = o.Int()
            is.Mb = o.Int()
        case strings.EqualFold(o.Name(), "ma"):
            is.Ma = o.Int()
        case strings.EqualFold(o.Name(), "mb"):
            is.Mb = o.Int()
        case strings.EqualFold(o.Name(), "k"):
            is.K = o.Int()
        case strings.EqualFold(o.Name(), "kl"):
            is.Kl = o.Int()
        case strings.EqualFold(o.Name(), "ku"):
            is.Ku = o.Int()
        case strings.EqualFold(o.Name(), "nrhs"):
            is.Nrhs = o.Int()
        }
    }
    return is
}

func PrintIndexes(p *IndexOpts) {
    // these used in BLAS/LAPACK
    fmt.Printf("N=%d, Nx=%d, Ny=%d\n", p.N, p.Nx, p.Ny)
    fmt.Printf("M=%d, Ma=%d, Mb=%d\n", p.M, p.Ma, p.Mb)
    fmt.Printf("LDa=%d, LDb=%d, LDc=%d\n", p.LDa, p.LDb, p.LDc)
    fmt.Printf("IncX=%d, IncY=%d\n", p.IncX, p.IncY)
    fmt.Printf("Ox=%d, Oy=%d, Oa=%d, Ob=%d, Oc=%d\n",
        p.OffsetX, p.OffsetY, p.OffsetA, p.OffsetB, p.OffsetC)
    fmt.Printf("K=%d, Ku=%d, Kl=%d\n", p.K, p.Ku, p.Kl)
    // these used in LAPACK
    fmt.Printf("NRHS=%d\n", p.Nrhs)
    fmt.Printf("Od=%d, Odl=%d, Odu=%d\n", p.OffsetD, p.OffsetDL, p.OffsetDU)
    fmt.Printf("LDw=%d, LDz=%d\n", p.LDw, p.LDz)
    fmt.Printf("Ow=%d, Oz=%d\n", p.OffsetW, p.OffsetZ)

}

// Local Variables:
// tab-width: 4
// End:
