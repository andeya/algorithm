// Copyright 2012 The geom Authors. All rights reserved.
// Use of this source code is governed by a license that
// can be found in the LICENSE file.

package geom

type Triangle struct {
	A, B, C Coord
}

func (me *Triangle) Bounds() (bounds *Rect) {
	bounds = &Rect{me.A, me.A}
	bounds.ExpandToContainCoord(me.B)
	bounds.ExpandToContainCoord(me.C)
	return
}

func (me *Triangle) Equals(oi interface{}) bool {
	ot, ok := oi.(*Triangle)
	if !ok {
		return false
	}
	if me.A.EqualsCoord(ot.A) {
		if me.B.EqualsCoord(ot.B) {
			return me.C.EqualsCoord(ot.C)
		}
		if me.B.EqualsCoord(ot.C) {
			return me.C.EqualsCoord(ot.B)
		}
	}
	if me.A.EqualsCoord(ot.B) {
		if me.B.EqualsCoord(ot.A) {
			return me.C.EqualsCoord(ot.C)
		}
		if me.B.EqualsCoord(ot.C) {
			return me.C.EqualsCoord(ot.A)
		}
	}
	if me.A.EqualsCoord(ot.C) {
		if me.B.EqualsCoord(ot.B) {
			return me.C.EqualsCoord(ot.A)
		}
		if me.B.EqualsCoord(ot.A) {
			return me.C.EqualsCoord(ot.B)
		}
	}
	return false
}

func (me *Triangle) Vertices() (vertices []Coord) {
	vertices = []Coord{me.A, me.B, me.C}
	return
}

func (me *Triangle) ContainsCoord(p Coord) bool {
	leftA := CrossProduct(me.B.Minus(me.A), p.Minus(me.A)) < 0
	leftB := CrossProduct(me.C.Minus(me.B), p.Minus(me.B)) < 0
	leftC := CrossProduct(me.A.Minus(me.C), p.Minus(me.C)) < 0
	return leftA == leftB && leftA == leftC
}

func (me *Triangle) HasVertex(v Coord) bool {
	return v.EqualsCoord(me.A) || v.EqualsCoord(me.B) || v.EqualsCoord(me.C)
}
