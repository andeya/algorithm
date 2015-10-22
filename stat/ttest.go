package stat

/* ttest.go
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

// runs a t-test between two datasets representing independent
// samples. Tests to see if the difference between means of the
// samples is different from zero.
func TTest(data1, data2 Interface) float64 {

	n1 := float64(data1.Len())
	n2 := float64(data2.Len())

	mean1 := Mean(data1)
	mean2 := Mean(data2)

	pv := PVariance(data1, data2)

	t := (mean1 - mean2) / (math.Sqrt(pv * ((1.0 / n1) + (1.0 / n2))))

	return t
}
