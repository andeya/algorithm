// Copyright 2012 The geom Authors. All rights reserved.
// Use of this source code is governed by a license that
// can be found in the LICENSE file.

package geom

import (
	"fmt"
	"math"
)

type Rect struct {
	Min, Max Coord
}

// this rect contains nothing
func NilRect() (r Rect) {
	r.Min.X = math.Inf(1)
	r.Min.Y = math.Inf(1)
	r.Max.X = math.Inf(-1)
	r.Max.Y = math.Inf(-1)
	return
}

func (r Rect) Width() float64 {
	return r.Max.X - r.Min.X
}

func (r Rect) Height() float64 {
	return r.Max.Y - r.Min.Y
}

func (r Rect) Size() (w, h float64) {
	return r.Max.X - r.Min.X, r.Max.Y - r.Min.Y
}

func (r Rect) Center() (center Coord) {
	center.X = 0.5 * (r.Min.X + r.Max.X)
	center.Y = 0.5 * (r.Min.Y + r.Max.Y)
	return
}

func (r Rect) ContainsCoord(p Coord) bool {
	return r.Min.QuadPP(p) && r.Max.QuadMM(p)
}

func (r Rect) ContainsRect(o Rect) bool {
	return r.ContainsCoord(o.Min) && r.ContainsCoord(o.Max)
}

func (r *Rect) Translate(offset Coord) {
	r.Min = r.Min.Plus(offset)
	r.Max = r.Max.Plus(offset)
}

func (r *Rect) Scale(xf, yf float64) {
	r.Min.Scale(xf, yf)
	r.Max.Scale(xf, yf)
	if xf < 0 {
		r.Min.X, r.Max.X = r.Max.X, r.Min.X
	}
	if yf < 0 {
		r.Min.Y, r.Max.Y = r.Max.Y, r.Min.Y
	}
}

func (r *Rect) ExpandToContain(ch <-chan Coord) {
	for p := range ch {
		r.ExpandToContainCoord(p)
	}
}

func (r *Rect) ExpandToContainCoord(p Coord) {
	r.Min.X = minf(r.Min.X, p.X)
	r.Min.Y = minf(r.Min.Y, p.Y)
	r.Max.X = maxf(r.Max.X, p.X)
	r.Max.Y = maxf(r.Max.Y, p.Y)
}

func (r *Rect) ExpandToContainRect(q Rect) {
	r.ExpandToContainCoord(q.Min)
	r.ExpandToContainCoord(q.Max)
}

func (r Rect) Bounds() (bounds Rect) {
	bounds = r
	return
}

func (r Rect) Equals(oi interface{}) bool {
	or, ok := oi.(Rect)
	return ok && RectsEqual(r, or)
}

func RectsIntersect(r1, r2 Rect) bool {
	ov := func(min1, max1, min2, max2 float64) (overlap bool) {
		if min1 <= min2 && max1 >= min2 {
			return true
		}
		if min1 <= max2 && max1 >= max2 {
			return true
		}
		if min2 <= min1 && max2 >= min1 {
			return true
		}
		if min2 <= max1 && max2 >= max1 {
			return true
		}
		return false
	}
	dbg("RI(%v, %v)", r1, r2)
	xoverlap := ov(r1.Min.X, r1.Max.X, r2.Min.X, r2.Max.X)
	yoverlap := ov(r1.Min.Y, r1.Max.Y, r2.Min.Y, r2.Max.Y)
	dbg("%v %v", xoverlap, yoverlap)
	return xoverlap && yoverlap
}

func RectsIntersectStrict(r1, r2 Rect) bool {
	ov := func(min1, max1, min2, max2 float64) (overlap bool) {
		if min1 < min2 && max1 > min2 {
			return true
		}
		if min1 < max2 && max1 > max2 {
			return true
		}
		if min2 < min1 && max2 > min1 {
			return true
		}
		if min2 < max1 && max2 > max1 {
			return true
		}
		return false
	}
	dbg("RI(%v, %v)", r1, r2)
	xoverlap := ov(r1.Min.X, r1.Max.X, r2.Min.X, r2.Max.X)
	yoverlap := ov(r1.Min.Y, r1.Max.Y, r2.Min.Y, r2.Max.Y)
	dbg("%v %v", xoverlap, yoverlap)
	return xoverlap && yoverlap
}

func RectsIntersection(r1, r2 Rect) (ri Rect) {
	ri.Min.X = math.Max(r1.Min.X, r2.Min.X)
	ri.Min.Y = math.Max(r1.Min.Y, r2.Min.Y)
	ri.Max.X = math.Min(r1.Max.X, r2.Max.X)
	ri.Max.Y = math.Min(r1.Max.Y, r2.Max.Y)
	return
}

func RectsEqual(r1, r2 Rect) bool {
	if !r1.Min.EqualsCoord(r2.Min) {
		return false
	}
	if !r1.Max.EqualsCoord(r2.Max) {
		return false
	}
	return true
}

func (r Rect) String() string {
	return fmt.Sprintf("{%v %v}", r.Min, r.Max)
}
