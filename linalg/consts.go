// Copyright (c) Harri Rautila, 2012

// This file is part of github.com/hrautila/linalg package. 
// It is free software, distributed under the terms of GNU Lesser General Public 
// License Version 3, or any later version. See the COPYING tile included in this archive.

package linalg

import (
    "errors"
    "fmt"
    "strings"
)

// BLAS/LAPACK matrix parameter constants.
const (
    // BLAS/LAPACK parameters. Chosen values match corresponding
    // parameters in CBLAS implementation.
    RowMajor    = 101 // Atlas row-major
    ColumnMajor = 102 // Atlas column major
    PNoTrans    = 111 // 'N'
    PTrans      = 112 // 'T'
    PConjTrans  = 113 // 'C'
    PUpper      = 121 // 'U'
    PLower      = 122 // 'L'
    PNonUnit    = 131 // 'N'
    PUnit       = 132 // 'U'
    PDiag       = 133 // 'D'
    PLeft       = 141 // 'L'
    PRight      = 142 // 'R'
    // These for LAPACK only
    PJobNo      = 151 // 'N'
    PJobValue   = 152 // 'V'
    PJobAll     = 153 // 'A'
    PJobS       = 154 // 'S'
    PJobO       = 155 // 'O'
    PRangeAll   = 161 // 'A'
    PRangeValue = 162 // 'V'
    PRangeInt   = 163 // 'I'
)

// Structure for BLAS/LAPACK function parameters.
type Parameters struct {
    Trans, TransA, TransB int
    Uplo                  int
    Diag                  int
    Side                  int
    Jobz                  int
    Jobu                  int
    Jobvt                 int
    Range                 int
}

func GetParam(name string, params ...Option) (val int) {
    val = -1
    for _, o := range params {
        if strings.EqualFold(o.Name(), name) {
            val = o.Int()
            return
        }
    }
    return
}

// Parse options and return parameter structure with option fields
// set to given or sensible defaults.
func GetParameters(params ...Option) (p *Parameters, err error) {
    err = nil
    p = &Parameters{
        PNoTrans,  // Trans
        PNoTrans,  // TransA
        PNoTrans,  // TransB
        PLower,    // Uplo
        PNonUnit,  // Diag
        PLeft,     // Side
        PJobNo,    // Jobz
        PJobNo,    // Jobu
        PJobNo,    // Jobvt
        PRangeAll} // Range

Loop:
    for _, o := range params {
        if _, ok := o.(*IOpt); !ok {
            continue Loop
        }
        pval := o.Int()
        switch {
        case strings.EqualFold(o.Name(), "trans"):
            if pval == PNoTrans || pval == PTrans || pval == PConjTrans {
                p.Trans = pval
                p.TransA = p.Trans
                p.TransB = p.Trans
            } else {
                err = errors.New("Illegal value for Transpose parameter")
                break Loop
            }
        case strings.EqualFold(o.Name(), "transa"):
            if pval == PNoTrans || pval == PTrans || pval == PConjTrans {
                p.TransA = pval
            } else {
                err = errors.New("Illegal value for Transpose parameter")
                break Loop
            }
        case strings.EqualFold(o.Name(), "transb"):
            if pval == PNoTrans || pval == PTrans || pval == PConjTrans {
                p.TransB = pval
            } else {
                err = errors.New("Illegal value for Transpose parameter")
                break Loop
            }
        case strings.EqualFold(o.Name(), "uplo"):
            if pval == PUpper || pval == PLower {
                p.Uplo = pval
            } else {
                err = errors.New("Illegal value for UpLo parameter")
                break Loop
            }
        case strings.EqualFold(o.Name(), "diag"):
            if pval == PNonUnit || pval == PUnit || pval == PDiag {
                p.Diag = pval
            } else {
                err = errors.New("Illegal value for Diag parameter")
                break Loop
            }
        case strings.EqualFold(o.Name(), "side"):
            if pval == PLeft || pval == PRight {
                p.Side = pval
            } else {
                err = errors.New("Illegal value for Side parameter")
                break Loop
            }
        // Lapack parameters
        case strings.EqualFold(o.Name(), "jobz"):
            if pval == PJobNo || pval == PJobValue {
                p.Jobz = pval
            } else {
                err = errors.New("Illegal value for Jobz parameter")
                break Loop
            }
        case strings.EqualFold(o.Name(), "jobu"):
            if pval == PJobNo || pval == PJobAll || pval == PJobS || pval == PJobO {
                p.Jobu = pval
            } else {
                err = errors.New("Illegal value for Jobu parameter")
                break Loop
            }
        case strings.EqualFold(o.Name(), "jobvt"):
            if pval == PJobNo || pval == PJobAll || pval == PJobS || pval == PJobO {
                p.Jobvt = pval
            } else {
                err = errors.New("Illegal value for Jobu parameter")
                break Loop
            }
        case strings.EqualFold(o.Name(), "range"):
            if pval == PRangeAll || pval == PRangeValue || pval == PRangeInt {
                p.Range = pval
            } else {
                err = errors.New("Illegal value for Range parameter")
                break Loop
            }
        }
    }
    return
}

// Matrix parameter option variables.
var (
    // trans: No Transpose 'N'
    OptNoTrans  = &IOpt{"trans", PNoTrans}
    OptNoTransA = &IOpt{"transA", PNoTrans}
    OptNoTransB = &IOpt{"transB", PNoTrans}
    // trans: Transpose 'T'
    OptTrans  = &IOpt{"trans", PTrans}
    OptTransA = &IOpt{"transA", PTrans}
    OptTransB = &IOpt{"transB", PTrans}
    // trans: Conjugate Transpose 'C'
    OptConjTrans  = &IOpt{"trans", PConjTrans}
    OptConjTransA = &IOpt{"transA", PConjTrans}
    OptConjTransB = &IOpt{"transB", PConjTrans}
    // uplo: Upper Triangular 'U'
    OptUpper = &IOpt{"uplo", PUpper}
    // uplo: Lower Triangular 'L'
    OptLower = &IOpt{"uplo", PLower}
    // side parameter 'L', 'R'
    OptLeft  = &IOpt{"side", PLeft}
    OptRight = &IOpt{"side", PRight}
    // diag parameter 'U'
    OptUnit = &IOpt{"diag", PUnit}
    // diag parameter 'N'
    OptNonUnit = &IOpt{"diag", PNonUnit}
    OptDiag    = &IOpt{"diag", PDiag}
    // Lapack jobz  'N'
    OptJobZNo = &IOpt{"jobz", PJobNo}
    // Lapack jobz   'V'
    OptJobZValue = &IOpt{"jobz", PJobValue}
    // Lapack jobu 'N'
    OptJobuNo = &IOpt{"jobu", PJobNo}
    // Lapack jobu 'A'
    OptJobuAll = &IOpt{"jobu", PJobAll}
    // Lapack jobu 'S'
    OptJobuS = &IOpt{"jobu", PJobS}
    // Lapack jobu 'O'
    OptJobuO = &IOpt{"jobu", PJobO}
    // Lapack jobvt 'N',
    OptJobvtNo = &IOpt{"jobvt", PJobNo}
    // Lapack jobvt 'A',
    OptJobvtAll = &IOpt{"jobvt", PJobAll}
    // Lapack jobvt 'S',
    OptJobvtS = &IOpt{"jobvt", PJobS}
    // Lapack jobvt 'O',
    OptJobvtO = &IOpt{"jobvt", PJobO}
    // Lapack range 'A'
    OptRangeAll = &IOpt{"range", PRangeAll}
    // Lapack range 'V'
    OptRangeValue = &IOpt{"range", PRangeValue}
    // Lapack range 'I'
    OptRangeInt = &IOpt{"range", PRangeInt}
)

var paramString map[int]string = map[int]string{
    PNoTrans:    "N",
    PTrans:      "T",
    PConjTrans:  "C",
    PUpper:      "U",
    PLower:      "L",
    PLeft:       "L",
    PRight:      "R",
    PUnit:       "U",
    PNonUnit:    "N",
    PJobNo:      "N",
    PJobValue:   "V",
    PJobAll:     "A",
    PJobS:       "S",
    PJobO:       "O",
    PRangeAll:   "A",
    PRangeValue: "V",
    PRangeInt:   "I",
}

// Map parameter value to name string that can be used when calling Fortran
// library functions.
func ParamString(p int) string {
    v, ok := paramString[p]
    if ok {
        return v
    }
    return ""
}

// Print parameter structure.
func PrintParameters(p *Parameters) {
    fmt.Printf("trans : %d [%s]\n", p.Trans, ParamString(p.Trans))
    fmt.Printf("transA: %d [%s]\n", p.TransA, ParamString(p.TransA))
    fmt.Printf("transB: %d [%s]\n", p.TransB, ParamString(p.TransB))
    fmt.Printf("Uplo  : %d [%s]\n", p.Uplo, ParamString(p.Uplo))
    fmt.Printf("Diag  : %d [%s]\n", p.Diag, ParamString(p.Diag))
    fmt.Printf("Side  : %d [%s]\n", p.Side, ParamString(p.Side))
    fmt.Printf("Jobz  : %d [%s]\n", p.Jobz, ParamString(p.Jobz))
    fmt.Printf("Range : %d [%s]\n", p.Range, ParamString(p.Range))
}

// Local Variables:
// tab-width: 4
// End:
