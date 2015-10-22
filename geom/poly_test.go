// Copyright 2012 The geom Authors. All rights reserved.
// Use of this source code is governed by a license that
// can be found in the LICENSE file.

package geom

import (
	"fmt"
	"testing"
)

func debugln(args ...interface{}) {
	if false {
		fmt.Println(args...)
	}
}
func debugf(format string, args ...interface{}) {
	if false {
		fmt.Printf(format, args...)
	}
}

func TestInsert(t *testing.T) {
	p := &Polygon{}
	p.AddVertex(Coord{0, 0})
	p.AddVertex(Coord{0, 1})
	p.AddVertex(Coord{1, 1})
	p.AddVertex(Coord{1, 0})
	p.InsertVertexAfter(Coord{5, 0}, 2)

	p2 := &Polygon{}
	p2.AddVertex(Coord{0, 0})
	p2.AddVertex(Coord{0, 1})
	p2.AddVertex(Coord{5, 0})
	p2.AddVertex(Coord{1, 1})
	p2.AddVertex(Coord{1, 0})

	if !p.Equals(p2) {
		t.Fail()
	}
}

func TestPolyTriangularize(t *testing.T) {
	poly := new(Polygon)
	poly.AddVertex(Coord{0, 0})
	poly.AddVertex(Coord{0, 1})
	poly.AddVertex(Coord{1, 1})
	poly.AddVertex(Coord{1, 0})
	tris, ok := poly.Triangles()
	if ok {
		debugln()
		for _, tri := range tris {
			debugf("triangle: %v\n", tri)
		}
	} else {
		debugf("No triangles for %v\n", poly)
	}

	poly = new(Polygon)
	poly.AddVertex(Coord{0, 0})
	poly.AddVertex(Coord{1, 1})
	poly.AddVertex(Coord{2, 0})
	poly.AddVertex(Coord{2, 3})
	poly.AddVertex(Coord{1, 2})
	poly.AddVertex(Coord{0, 3})
	tris, ok = poly.Triangles()
	if ok {
		debugln()
		for _, tri := range tris {
			debugf("triangle: %v\n", tri)
		}
	} else {
		debugf("No triangles for %v\n", poly)
	}

	poly = new(Polygon)
	poly.AddVertex(Coord{2, 1})
	poly.AddVertex(Coord{2, 2})
	poly.AddVertex(Coord{1, 2})
	poly.AddVertex(Coord{1, 3})
	tris, ok = poly.Triangles()
	if ok {
		debugln()
		for _, tri := range tris {
			debugf("triangle: %v\n", tri)
		}
	} else {
		debugf("No triangles for %v\n", poly)
	}
}

//{44 736} {44 848} {88 848} {88 1044} {44 1044} {44 1244} {68 1244} {68 1068} {112 1068} {112 824} {68 824} {68 736}
func TestPiece(t *testing.T) {
	vertices := []Coord{
		Coord{1, 1},
		Coord{1, 6},
		Coord{2, 6},
		Coord{2, 3},
		Coord{4, 3},
		Coord{4, 6},
		Coord{5, 6},
		Coord{5, 2},
		Coord{2, 2},
		Coord{2, 1},
	}

	poly := new(Polygon)
	for _, v := range vertices {
		poly.AddVertex(v)
	}
	tris, ok := poly.Triangles()
	if ok {
		debugln()
		for _, tri := range tris {
			debugf("triangle: %v\n", tri)
		}
	} else {
		debugf("No triangles for %v\n", poly)
	}

	vertices = []Coord{
		Coord{44, 736},
		Coord{44, 848},
		Coord{88, 848},
		Coord{88, 1044},
		Coord{44, 1044},
		Coord{44, 1244},
		Coord{68, 1244},
		Coord{68, 1068},
		Coord{112, 1068},
		Coord{112, 824},
		Coord{68, 824},
		Coord{68, 736},
	}

	poly = new(Polygon)
	for _, v := range vertices {
		poly.AddVertex(v)
	}
	tris, ok = poly.Triangles()
	if ok {
		debugln()
		for _, tri := range tris {
			debugf("triangle: %v\n", tri)
		}
	} else {
		debugf("No triangles for %v\n", poly)
	}
}
