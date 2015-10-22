// Copyright 2012 The geom Authors. All rights reserved.
// Use of this source code is governed by a license that
// can be found in the LICENSE file.

package geom

type Bounded interface {
	Bounds() (bounds Rect)
}

type Transformable interface {
	Translate(offset Coord)
	Rotate(rad float64)
	Scale(xfactor, yfactor float64)
}
