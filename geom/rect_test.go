// Copyright 2012 The geom Authors. All rights reserved.
// Use of this source code is governed by a license that
// can be found in the LICENSE file.

package geom

import (
	"testing"
)

func TestRectsIntersect(t *testing.T) {
	r1 := Rect{Coord{0, 0}, Coord{650, 650}}
	r2 := Rect{Coord{200, 500}, Coord{450, 750}}
	if !r1.Min.QuadPP(r2.Min) {
		t.Error("QuadPP")
	}
	if !r1.Max.QuadMM(r2.Min) {
		t.Error("QuadMM")
	}
	if !r1.ContainsCoord(r2.Min) {
		t.Error("contains")
	}
	if !RectsIntersect(r1, r2) {
		t.Error("intersect")
	}

	r1 = Rect{Coord{325, 325}, Coord{650, 650}}
	r2 = Rect{Coord{200, 500}, Coord{450, 750}}

	Debug = false

	if !r1.Min.QuadMP(r2.Min) {
		t.Error("QuadMP2")
	}
	if !r1.Max.QuadMM(r2.Min) {
		t.Error("QuadMM2")
	}
	if r1.ContainsCoord(r2.Min) {
		t.Error("contains2")
	}
	if !RectsIntersect(r1, r2) {
		t.Error("intersect2")
	}

}
