// Copyright 2012 The geom Authors. All rights reserved.
// Use of this source code is governed by a license that
// can be found in the LICENSE file.

package geom

func minf(x, y float64) (r float64) {
	r = x
	if r > y {
		r = y
	}
	return
}

func maxf(x, y float64) (r float64) {
	r = x
	if r < y {
		r = y
	}
	return
}
