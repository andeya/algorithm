// Copyright 2012 The geom Authors. All rights reserved.
// Use of this source code is governed by a license that
// can be found in the LICENSE file.

package qtree

import (
	"fmt"
	"github.com/henrylee2cn/algorithm/geom"
	"testing"
)

func TestInsertCollect(t *testing.T) {
	cfg := ConfigDefault()
	qt := New(cfg, geom.Rect{geom.Coord{0, 0}, geom.Coord{100, 100}})

	r := geom.Rect{geom.Coord{20, 20}, geom.Coord{40, 40}}
	qt.Insert(r)

	collection := make(map[Item]bool)
	qt.CollectIntersect(r, collection)

	fmt.Printf("%v\n", collection)
}
