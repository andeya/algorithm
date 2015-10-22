// Copyright 2012 The geom Authors. All rights reserved.
// Use of this source code is governed by a license that
// can be found in the LICENSE file.

package geom

/*
A line that goes through Intersection and has a normal Normal.
*/
type Line struct {
	Intersection Coord
	Normal       Coord
}

func LineIntersection(l1, l2 Line) (p Coord) {
	b1 := l1.Normal.Unit()
	b1.RotateRight()
	unum := (l2.Normal.X*l2.Intersection.X + l2.Normal.Y*l2.Intersection.Y)
	uden := (l2.Normal.X*b1.X + l2.Normal.Y*b1.Y)
	u := unum / uden
	p = l1.Intersection.Plus(b1.Times(u))
	return
}
