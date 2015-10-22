package stat

/* wabsdev.go
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Jim Davies, Brian Gough
 * Copyright (C) 2012, 2013 G.vd.Schoot
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

import (
	"math"
)

func WAbsdev(w, data Interface) float64 {
	wmean := WMean(w, data)
	return WAbsdevMean(w, data, wmean)
}

// WAbsdevMean calculates the weighted absolute deviation of a dataset
func WAbsdevMean(w, data Interface, wmean float64) (wabsdev float64) {

	var W float64
	Len := data.Len()

	// calculate the sum of the absolute deviations
	for i := 0; i < Len; i++ {
		wi := w.Get(i)

		if wi > 0 {
			delta := math.Abs(data.Get(i) - wmean)
			W += wi
			wabsdev += (delta - wabsdev) * (wi / W)
		}
	}

	return
}
