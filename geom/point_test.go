// Copyright 2012 The geom Authors. All rights reserved.
// Use of this source code is governed by a license that
// can be found in the LICENSE file.

package geom

import (
	"fmt"
	"testing"
)

func TestVertexAngle(t *testing.T) {
	A := Coord{0, 0}
	B := Coord{1, 0}
	C := Coord{1, 1}
	D := Coord{1, -1}
	r1 := VertexAngle(A, B, C)
	r2 := VertexAngle(A, B, D)
	if r1 != -r2 {

	}

	p1 := Coord{1, 2}
	p2 := Coord{0, 3}
	p3 := Coord{0, 0}
	p4 := Coord{1, 1}

	rp := VertexAngle(p1, p2, p3)
	rn := VertexAngle(p4, p2, p3)
	fmt.Println(rp, rn)
}

func TestVectorAngle(t *testing.T) {
	v1 := Coord{0, -1}
	v2 := Coord{-1, 0}
	v3 := Coord{1, -1}
	v4 := Coord{-1, 0}

	a12 := VectorAngle(v1, v2)
	a34 := VectorAngle(v3, v4)

	fmt.Println(a12, a34)
}
