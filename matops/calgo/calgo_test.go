// Copyright (c) Harri Rautila, 2012,2013

// This file is part of github.com/hrautila/matops package. It is free software,
// distributed under the terms of GNU Lesser General Public License Version 3, or
// any later version. See the COPYING tile included in this archive.

package calgo

import (
	"github.com/henrylee2cn/algorithm/linalg/blas"
	"github.com/henrylee2cn/algorithm/matrix"
	//"github.com/henrylee2cn/algorithm/linalg"
	"math"
	"math/rand"
	"testing"
	"time"
)

const M = 8
const N = 8
const P = 8

var data7 [][]float64 = [][]float64{
	[]float64{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
	[]float64{2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0},
	[]float64{3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0},
	[]float64{4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0},
	[]float64{5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0},
	[]float64{6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0},
	[]float64{7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0}}

var data6 [][]float64 = [][]float64{
	[]float64{1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
	[]float64{2.0, 2.0, 2.0, 2.0, 2.0, 2.0},
	[]float64{3.0, 3.0, 3.0, 3.0, 3.0, 3.0},
	[]float64{4.0, 4.0, 4.0, 4.0, 4.0, 4.0},
	[]float64{5.0, 5.0, 5.0, 5.0, 5.0, 5.0},
	[]float64{6.0, 6.0, 6.0, 6.0, 6.0, 6.0}}

var data5 [][]float64 = [][]float64{
	[]float64{1.0, 1.0, 1.0, 1.0, 1.0},
	[]float64{2.0, 2.0, 2.0, 2.0, 2.0},
	[]float64{3.0, 3.0, 3.0, 3.0, 3.0},
	[]float64{4.0, 4.0, 4.0, 4.0, 4.0},
	[]float64{5.0, 5.0, 5.0, 5.0, 5.0}}

var upper7 [][]float64 = [][]float64{
	[]float64{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
	[]float64{0.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0},
	[]float64{0.0, 0.0, 3.0, 3.0, 3.0, 3.0, 3.0},
	[]float64{0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 4.0},
	[]float64{0.0, 0.0, 0.0, 0.0, 5.0, 5.0, 5.0},
	[]float64{0.0, 0.0, 0.0, 0.0, 0.0, 6.0, 6.0},
	[]float64{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.0}}

var lower7 [][]float64 = [][]float64{
	[]float64{1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	[]float64{1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	[]float64{1.0, 2.0, 3.0, 0.0, 0.0, 0.0, 0.0},
	[]float64{1.0, 2.0, 3.0, 4.0, 0.0, 0.0, 0.0},
	[]float64{1.0, 2.0, 3.0, 4.0, 5.0, 0.0, 0.0},
	[]float64{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 0.0},
	[]float64{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0}}

var lower5 [][]float64 = [][]float64{
	[]float64{1.0, 0.0, 0.0, 0.0, 0.0},
	[]float64{1.0, 2.0, 0.0, 0.0, 0.0},
	[]float64{1.0, 2.0, 3.0, 0.0, 0.0},
	[]float64{1.0, 2.0, 3.0, 4.0, 0.0},
	[]float64{1.0, 2.0, 3.0, 4.0, 5.0}}

var upper3 [][]float64 = [][]float64{
	[]float64{1.0, 1.0, 1.0},
	[]float64{0.0, 2.0, 2.0},
	[]float64{0.0, 0.0, 3.0}}

var upper3_1 [][]float64 = [][]float64{
	[]float64{2.0, 2.0, 2.0},
	[]float64{0.0, 3.0, 3.0},
	[]float64{0.0, 0.0, 4.0}}

var lower3 [][]float64 = [][]float64{
	[]float64{1.0, 0.0, 0.0},
	[]float64{1.0, 2.0, 0.0},
	[]float64{1.0, 2.0, 3.0}}

var data7x2_iinc [][]float64 = [][]float64{
	[]float64{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0},
	[]float64{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0}}

var data7x2_dec [][]float64 = [][]float64{
	[]float64{7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0},
	[]float64{7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0}}

func fail(t *testing.T, isok bool) {
	if !isok {
		t.Errorf("** FAILURE HERE ***\n")
	}
}

func TestMakeData(t *testing.T) {
	blas.PanicOnError(true)
	rand.Seed(time.Now().UnixNano())
}

func isClose(a, b float64) bool {
	const RTOL = 1.0000000000000001e-05
	const ATOL = 1e-8
	df := math.Abs(a - b)
	if df > ATOL+RTOL*math.Abs(b) {
		return false
	}
	return true
}

func TestVectors(t *testing.T) {
	a := matrix.FloatWithValue(1, 5, 1.0)
	b := matrix.FloatWithValue(1, 5, 2.0)
	c := matrix.FloatWithValue(5, 1, 2.0)
	ar := a.FloatArray()
	br := b.FloatArray()
	cr := c.FloatArray()
	v := DDot(ar, br, 1.0, a.LeadingIndex(), b.LeadingIndex(), a.NumElements())
	t.Logf("a*b = %.1f\n", v)
	v = DDot(cr, br, 1.0, 1, b.LeadingIndex(), c.NumElements())
	t.Logf("c*b = %.1f\n", v)
	v = DNorm2(br, 1, b.NumElements())
	t.Logf("norm2(b) = %.1f\n", v)

	b0 := matrix.FloatNew(1, 5, []float64{-1.0, 0.0, 2.0, -3.0, 5.0})
	br = b0.FloatArray()
	ix := DIAMax(br, 1, b0.NumElements())
	t.Logf("iamax = %d: %.1f %v\n", ix, b0.GetIndex(ix), b0)

	b0 = matrix.FloatNew(1, 5, []float64{-8.0, 0.0, 2.0, -3.0, 5.0})
	br = b0.FloatArray()
	ix = DIAMax(br, 1, b0.NumElements())
	t.Logf("iamax = %d: %.1f %v\n", ix, b0.GetIndex(ix), b0)

	b0 = matrix.FloatNew(1, 5, []float64{-8.0, 0.0, 2.0, -9.0, 5.0})
	br = b0.FloatArray()
	ix = DIAMax(br, 1, b0.NumElements())
	t.Logf("iamax = %d: %.1f %v\n", ix, b0.GetIndex(ix), b0)
}

// Local Variables:
// tab-width: 4
// indent-tabs-mode: nil
// End:
