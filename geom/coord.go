// Copyright 2012 The geom Authors. All rights reserved.
// Use of this source code is governed by a license that
// can be found in the LICENSE file.

package geom

import (
	"math"
)

type Coord struct {
	X, Y float64
}

func (p *Coord) Hashcode() (hash uint64) {
	x, y := uint64(p.X), uint64(p.Y)
	hash = x + y
	return
}

func (p *Coord) Equals(oi interface{}) (equals bool) {
	o, equals := oi.(*Coord)
	if !equals {
		var op Coord
		op, equals = oi.(Coord)
		equals = equals && p.EqualsCoord(op)
		return
	}
	equals = p.EqualsCoord(*o)
	return
}

func (p *Coord) Translate(offset Coord) {
	*p = p.Plus(offset)
}

func (p *Coord) Rotate(rad float64) {
	p.X = p.X*math.Cos(rad) - p.Y*math.Sin(rad)
	p.Y = p.X*math.Sin(rad) + p.Y*math.Cos(rad)
}

func (p *Coord) RotateLeft() {
	p.X, p.Y = -p.Y, p.X
}

func (p *Coord) RotateRight() {
	p.X, p.Y = p.Y, -p.X
}

func (p Coord) Unit() (u Coord) {
	m := p.Magnitude()
	u.X = p.X / m
	u.Y = p.Y / m
	return
}

func (p *Coord) Scale(xfactor, yfactor float64) {
	p.X *= xfactor
	p.Y *= yfactor
}

func (p Coord) EqualsCoord(q Coord) bool {
	return p.X == q.X && p.Y == q.Y
}

func (p Coord) DistanceFrom(q Coord) (d float64) {
	return p.Minus(q).Magnitude()
}

func (p Coord) DistanceFromSquared(q Coord) (ds float64) {
	return p.Minus(q).MagnitudeSquared()
}

func (p Coord) Magnitude() (m float64) {
	m = math.Sqrt(p.MagnitudeSquared())
	return
}

func (p Coord) MagnitudeSquared() (ms float64) {
	ms = p.X*p.X + p.Y*p.Y
	return
}

func (p Coord) Minus(q Coord) (r Coord) {
	r.X = p.X - q.X
	r.Y = p.Y - q.Y
	return
}

func (p Coord) Plus(q Coord) (r Coord) {
	r.X = p.X + q.X
	r.Y = p.Y + q.Y
	return
}

func (p Coord) Times(s float64) (r Coord) {
	r.X = p.X * s
	r.Y = p.Y * s
	return
}

func (p Coord) QuadPP(q Coord) bool {
	return q.X >= p.X && q.Y >= p.Y
}

func (p Coord) QuadPM(q Coord) bool {
	return q.X >= p.X && q.Y <= p.Y
}

func (p Coord) QuadMP(q Coord) bool {
	return q.X <= p.X && q.Y >= p.Y
}

func (p Coord) QuadMM(q Coord) bool {
	return q.X <= p.X && q.Y <= p.Y
}

func DotProduct(p, q Coord) (r float64) {
	r = p.X*q.X + p.Y*q.Y
	return
}

func CrossProduct(p, q Coord) (z float64) {
	z = p.X*q.Y - p.Y*q.X
	return
}

func VectorAngle(X, Y Coord) (r float64) {
	XdotY := DotProduct(X, Y)
	mXmY := X.Magnitude() * Y.Magnitude()
	r = math.Acos(XdotY / mXmY)
	z := CrossProduct(X, Y)
	if z < 0 {
		r *= -1
	}
	return
}

func VertexAngle(A, B, C Coord) (r float64) {
	X := A.Minus(B)
	Y := C.Minus(B)
	r = VectorAngle(X, Y)
	if r < 0 {
		r *= -1
	}
	return
}

func CoordChan(points []Coord) (ch <-chan Coord) {
	tch := make(chan Coord, len(points))
	go func(points []Coord, ch chan<- Coord) {
		for _, p := range points {
			ch <- p
		}
		close(ch)
	}(points, tch)
	ch = tch
	return
}
