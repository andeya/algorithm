// Copyright 2012 The geom Authors. All rights reserved.
// Use of this source code is governed by a license that
// can be found in the LICENSE file.

package geom

import (
	"math"
)

type Polygon struct {
	Path
}

func wrapIndex(index, length int) (i int) {
	i = index % length
	if i < 0 {
		i = length + i
	}
	return
}

func (p *Polygon) Clone() (op *Polygon) {
	op = &Polygon{*p.Path.Clone()}
	return
}

func (p *Polygon) Equals(oi interface{}) bool {
	o, ok := oi.(*Polygon)
	if !ok {
		return false
	}
	return (&p.Path).Equals(&o.Path)
}

func (p *Polygon) Register(op *Polygon) (offset Coord, match bool) {
	offset, match = p.Path.Register(&op.Path)
	return
}

func (me *Polygon) Vertex(index int) (v Coord) {
	v = me.vertices[wrapIndex(index, len(me.vertices))]
	return
}

func (me *Polygon) Segment(index int) (s *Segment) {
	s = &Segment{me.Vertex(index), me.Vertex(index + 1)}
	return
}

func (me *Polygon) VertexAngle(index int) (r float64) {
	a := me.Vertex(index - 1)
	b := me.Vertex(index)
	c := me.Vertex(index + 1)
	r = VertexAngle(a, b, c)
	return
}

func (me *Polygon) WindingOrder() (winding float64) {
	for i := 0; i < len(me.vertices); i++ {
		winding += me.VertexAngle(i)
	}
	return
}

func (me *Polygon) ContainsCoord(p Coord) bool {
	fakeSegment := &Segment{p, Coord{p.X, p.Y + 1}}

	above := 0
	for i := 0; i < me.Length(); i++ {
		s := me.Segment(i)
		uh, uv := s.IntersectParameters(fakeSegment)
		if uh < 0 || uh >= 1 {
			continue
		}
		if uv > 0 {
			above++
		}
	}
	return above%2 == 1
}

//bisect a polygon by joining vertices i and j
func (me *Polygon) Bisect(i, j int) (p1, p2 *Polygon) {
	i = wrapIndex(i, len(me.vertices))
	j = wrapIndex(j, len(me.vertices))

	//build the first one, starting at i and ending at j
	p1 = &Polygon{}
	for c := i; c != wrapIndex(j+1, len(me.vertices)); c = wrapIndex(c+1, len(me.vertices)) {
		p1.AddVertex(me.Vertex(c))
	}

	//build the second one, starting at j and ending at i
	p2 = &Polygon{}
	for c := j; c != wrapIndex(i+1, len(me.vertices)); c = wrapIndex(c+1, len(me.vertices)) {
		p2.AddVertex(me.Vertex(c))
	}

	return
}

func (me *Polygon) Error(other *Polygon) (offset Coord, error float64) {
	return me.Path.Error(&other.Path)
}

func (me *Polygon) Triangles() (tris []Triangle, ok bool) {
	dbg("%v.Triangles()", me)

	if me.Length() == 3 {
		dbg("already a triangle")
		tris = []Triangle{Triangle{me.Vertex(0), me.Vertex(1), me.Vertex(2)}}
		ok = true
		return
	}

	for i := 0; i < me.Length(); i++ {
		iv := me.Vertex(i)
	v2:
		for j := i + 2; j != wrapIndex(i-1, me.Length()); j = wrapIndex(j+1, me.Length()) {
			jv := me.Vertex(j)
			bisectingSegment := &Segment{iv, jv}
			dbg("bisectingSegment(%d, %d) = %v", i, j, bisectingSegment)

			//first check to see that it doesn't intersect any other segments
			for si := 0; si < me.Length(); si++ {
				s := me.Segment(si)
				u1, u2 := s.IntersectParameters(bisectingSegment)
				if math.IsNaN(u1) || math.IsNaN(u2) || (u1 > 0 && u1 < 1 && u2 > 0 && u2 < 1) {
					dbg(" Segment(%d, %d) %v\n%f %f", si, si+1, s, u1, u2)
					continue v2
				} else {
					dbg(" doesn't intersect %v: %f %f", s, u1, u2)
				}
			}

			//second check to see that it is in the interior of the polygon
			midCoord := bisectingSegment.Extrapolate(0.5)
			if !me.ContainsCoord(midCoord) {
				dbg(" poly contains %v", midCoord)
				continue v2
			}

			dbg(" Segment %v is good", bisectingSegment)

			p1, p2 := me.Bisect(i, j)
			t1, ok1 := p1.Triangles()
			t2, ok2 := p2.Triangles()
			tris = append(t1, t2...)
			ok = ok1 && ok2
			return
		}
	}

	dbg("failed with %v", me)
	//panic("couldn't find any valid bisecting segment")

	return
}
