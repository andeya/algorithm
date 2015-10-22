henrylee2cn/algorithm// Copyright (c) Harri Rautila, 2012

// This file is part of github.com/hrautila/cvx package. 
// It is free software, distributed under the terms of GNU Lesser General Public 
// License Version 3, or any later version. See the COPYING tile included in this archive.

package cvx

import (
    "errors"
    "fmt"
    "github.com/henrylee2cn/algorithm/cvx/sets"
    la "github.com/henrylee2cn/algorithm/linalg"
    "github.com/henrylee2cn/algorithm/linalg/blas"
    "github.com/henrylee2cn/algorithm/matrix"
    "math"
)

func checkArgK(K []int) (err error) {
    err = nil
    if len(K) == 0 {
        err = errors.New("'K' must be a non-empty list of positive integers")
        return
    }
    for _, k := range K {
        if k <= 0 {
            err = errors.New("'K' must be a non-empty list of positive integers")
            return
        }
    }
    return
}

type gpConvexProg struct {
    mnl  int
    l    int
    n    int
    ind  [][2]int
    F    *matrix.FloatMatrix
    g    *matrix.FloatMatrix
    maxK int
}

func createGpProg(K []int, F, g *matrix.FloatMatrix) *gpConvexProg {
    gp := &gpConvexProg{mnl: len(K) - 1, l: F.Rows(), n: F.Cols()}
    gp.ind = make([][2]int, len(K))
    s := 0
    for i := 0; i < len(K); i++ {
        gp.ind[i][0] = s
        gp.ind[i][1] = s + K[i]
        s += K[i]
    }
    gp.F = F
    gp.g = g
    gp.maxK = maxdim(K)
    return gp
}

func (gp *gpConvexProg) F0() (mnl int, x *matrix.FloatMatrix, err error) {
    return gp.mnl, matrix.FloatZeros(gp.n, 1), nil
}

func (gp *gpConvexProg) F1(x *matrix.FloatMatrix) (f, Df *matrix.FloatMatrix, err error) {
    f = nil
    Df = nil
    err = nil
    f = matrix.FloatZeros(gp.mnl+1, 1)
    Df = matrix.FloatZeros(gp.mnl+1, gp.n)
    y := gp.g.Copy()
    blas.GemvFloat(gp.F, x, y, 1.0, 1.0)

    for i, s := range gp.ind {
        start := s[0]
        stop := s[1]
        // yi := exp(yi) = exp(Fi*x+gi)
        ymax := maxvec(y.FloatArray()[start:stop])
        // ynew = exp(y[start:stop] - ymax)
        ynew := matrix.Exp(matrix.FloatVector(y.FloatArray()[start:stop]).Add(-ymax))
        y.SetIndexesFromArray(ynew.FloatArray(), matrix.Indexes(start, stop)...)

        // fi = log sum yi = log sum exp(Fi*x+gi)
        ysum := blas.AsumFloat(y, &la.IOpt{"n", stop - start}, &la.IOpt{"offset", start})
        f.SetIndex(i, ymax+math.Log(ysum))

        blas.ScalFloat(y, 1.0/ysum, &la.IOpt{"n", stop - start}, &la.IOpt{"offset", start})
        blas.GemvFloat(gp.F, y, Df, 1.0, 0.0, la.OptTrans, &la.IOpt{"m", stop - start},
            &la.IOpt{"incy", gp.mnl + 1}, &la.IOpt{"offseta", start},
            &la.IOpt{"offsetx", start}, &la.IOpt{"offsety", i})
    }
    return
}

func (gp *gpConvexProg) F2(x, z *matrix.FloatMatrix) (f, Df, H *matrix.FloatMatrix, err error) {

    err = nil
    f = matrix.FloatZeros(gp.mnl+1, 1)
    Df = matrix.FloatZeros(gp.mnl+1, gp.n)
    H = matrix.FloatZeros(gp.n, gp.n)
    y := gp.g.Copy()
    Fsc := matrix.FloatZeros(gp.maxK, gp.n)
    blas.GemvFloat(gp.F, x, y, 1.0, 1.0)
    //fmt.Printf("y=\n%v\n", y.ToString("%.3f"))

    for i, s := range gp.ind {
        start := s[0]
        stop := s[1]

        // yi := exp(yi) = exp(Fi*x+gi)
        ymax := maxvec(y.FloatArray()[start:stop])
        ynew := matrix.Exp(matrix.FloatVector(y.FloatArray()[start:stop]).Add(-ymax))
        y.SetIndexesFromArray(ynew.FloatArray(), matrix.Indexes(start, stop)...)

        // fi = log sum yi = log sum exp(Fi*x+gi)
        ysum := blas.AsumFloat(y, &la.IOpt{"n", stop - start}, &la.IOpt{"offset", start})

        f.SetIndex(i, ymax+math.Log(ysum))
        blas.ScalFloat(y, 1.0/ysum, &la.IOpt{"n", stop - start}, &la.IOpt{"offset", start})
        blas.GemvFloat(gp.F, y, Df, 1.0, 0.0, la.OptTrans, &la.IOpt{"m", stop - start},
            &la.IOpt{"incy", gp.mnl + 1}, &la.IOpt{"offseta", start},
            &la.IOpt{"offsetx", start}, &la.IOpt{"offsety", i})

        Fsc.SetSubMatrix(0, 0, gp.F.GetSubMatrix(start, 0, stop-start))

        for k := start; k < stop; k++ {
            blas.AxpyFloat(Df, Fsc, -1.0, &la.IOpt{"n", gp.n},
                &la.IOpt{"incx", gp.mnl + 1}, &la.IOpt{"incy", Fsc.Rows()},
                &la.IOpt{"offsetx", i}, &la.IOpt{"offsety", k - start})
            blas.ScalFloat(Fsc, math.Sqrt(y.GetIndex(k)),
                &la.IOpt{"inc", Fsc.Rows()}, &la.IOpt{"offset", k - start})
        }
        // H += z[i]*Hi = z[i] *Fisc' * Fisc
        blas.SyrkFloat(Fsc, H, z.GetIndex(i), 1.0, la.OptTrans,
            &la.IOpt{"k", stop - start})
    }
    return
}

//
// Solves a geometric program
//
//   minimize    log sum exp (F0*x+g0)
//   subject to  log sum exp (Fi*x+gi) <= 0,  i=1,...,m
//               G*x <= h      
//               A*x = b
//
func Gp(K []int, F, g, G, h, A, b *matrix.FloatMatrix, solopts *SolverOptions) (sol *Solution, err error) {

    if err = checkArgK(K); err != nil {
        return
    }
    l := sumdim(K)

    if F == nil || F.Rows() != l {
        err = errors.New(fmt.Sprintf("'F' must matrix with %d rows", l))
        return
    }

    if g == nil || !g.SizeMatch(l, 1) {
        err = errors.New(fmt.Sprintf("'g' must matrix with size (%d,1)", l))
        return
    }
    n := F.Cols()

    if G == nil {
        G = matrix.FloatZeros(0, n)
    }
    if h == nil {
        h = matrix.FloatZeros(0, 1)
    }
    if G.Cols() != n {
        err = errors.New(fmt.Sprintf("'G' must matrix with size %d columns", n))
        return
    }
    ml := G.Rows()
    if h == nil || !h.SizeMatch(ml, 1) {
        err = errors.New(fmt.Sprintf("'h' must matrix with size (%d,1)", ml))
        return
    }

    if A == nil {
        A = matrix.FloatZeros(0, n)
    }
    if b == nil {
        b = matrix.FloatZeros(0, 1)
    }
    if A.Cols() != n {
        err = errors.New(fmt.Sprintf("'A' must matrix with size %d columns", n))
        return
    }
    p := A.Rows()
    if b == nil || !b.SizeMatch(p, 1) {
        err = errors.New(fmt.Sprintf("'b' must matrix with size (%d,1)", p))
        return
    }

    dims := sets.NewDimensionSet("l", "q", "s")
    dims.Set("l", []int{ml})
    gpProg := createGpProg(K, F, g)

    return Cp(gpProg, G, h, A, b, dims, solopts)
}

// Local Variables:
// tab-width: 4
// End:
