// Copyright 2012 The geom Authors. All rights reserved.
// Use of this source code is governed by a license that
// can be found in the LICENSE file.

package geom

type Segment struct {
	A, B Coord
}

func (s *Segment) IntersectParameters(t *Segment) (ps, pt float64) {
	x1, y1 := s.A.X, s.A.Y
	x2, y2 := s.B.X, s.B.Y
	x3, y3 := t.A.X, t.A.Y
	x4, y4 := t.B.X, t.B.Y

	x4x3 := x4 - x3
	y1y3 := y1 - y3
	y4y3 := y4 - y3
	x1x3 := x1 - x3
	x2x1 := x2 - x1
	y2y1 := y2 - y1

	ps = (x4x3*y1y3 - y4y3*x1x3) / (y4y3*x2x1 - x4x3*y2y1)
	pt = (x2x1*y1y3 - y2y1*x1x3) / (y4y3*x2x1 - x4x3*y2y1)

	return
}

func (s *Segment) Intersection(t *Segment) (p Coord, ok bool) {
	ps, pt := s.IntersectParameters(t)

	p.X = s.A.X + ps*(s.B.X-s.A.X)
	p.Y = s.A.Y + ps*(s.B.Y-s.A.Y)

	ok = ps >= 0 && ps <= 1 && pt >= 0 && pt <= 1

	return
}

func (s *Segment) Extrapolate(t float64) (p Coord) {
	p.X = s.A.X + t*(s.B.X-s.A.X)
	p.Y = s.A.Y + t*(s.B.Y-s.A.Y)
	return
}
