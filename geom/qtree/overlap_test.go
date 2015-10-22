// Copyright 2012 The geom Authors. All rights reserved.
// Use of this source code is governed by a license that
// can be found in the LICENSE file.

package qtree

/*
A number of tests to test qtree's ability to analyze the differences between
collections of shapes
*/

import (
	"github.com/henrylee2cn/algorithm/geom"
)

func col1() (items []Item) {
	items = append(items, Item(&geom.Rect{geom.Coord{0, 0}, geom.Coord{1, 1}}))
	items = append(items, Item(&geom.Rect{geom.Coord{0, 2}, geom.Coord{1, 3}}))
	items = append(items, Item(&geom.Rect{geom.Coord{2, 2}, geom.Coord{3, 3}}))
	items = append(items, Item(&geom.Rect{geom.Coord{2, 0}, geom.Coord{3, 1}}))
	return
}

func col2() (items []Item) {
	items = append(items, Item(&geom.Rect{geom.Coord{10, 0}, geom.Coord{11, 1}}))
	items = append(items, Item(&geom.Rect{geom.Coord{10, 2}, geom.Coord{11, 3}}))
	items = append(items, Item(&geom.Rect{geom.Coord{12, 2}, geom.Coord{13, 3}}))
	items = append(items, Item(&geom.Rect{geom.Coord{12, 0}, geom.Coord{13, 1}}))
	return
}

func col3() (items []Item) {
	items = append(items, Item(&geom.Rect{geom.Coord{0, 10}, geom.Coord{1, 11}}))
	items = append(items, Item(&geom.Rect{geom.Coord{0, 12}, geom.Coord{1, 13}}))
	items = append(items, Item(&geom.Rect{geom.Coord{2, 12}, geom.Coord{3, 13}}))
	items = append(items, Item(&geom.Rect{geom.Coord{2, 10}, geom.Coord{3, 11}}))
	return
}
