package stat

/* minmax.go
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

// Max finds the first largest member and the members position within the data 
func Max(data Interface) (max float64, max_index int) {
	Len := data.Len()
	max = data.Get(0)

	for i := 0; i < Len; i++ {
		xi := data.Get(i)

		if xi > max {
			max = xi
			max_index = i
		}

		if math.IsNaN(xi) {
			max = xi
			max_index = i
			return
		}
	}

	return
}

// Min finds the first smallest member and the members position within the data 
func Min(data Interface) (min float64, min_index int) {
	Len := data.Len()
	min = data.Get(0)

	for i := 0; i < Len; i++ {
		xi := data.Get(i)

		if xi < min {
			min = xi
			min_index = i
		}

		if math.IsNaN(xi) {
			min = xi
			min_index = i
			return
		}
	}

	return
}

// Minmax finds the first smallest and largest members and 
// the members positions within the data 
func Minmax(data Interface) (min float64, min_index int, max float64, max_index int) {
	Len := data.Len()
	min = data.Get(0)
	max = data.Get(0)

	for i := 0; i < Len; i++ {
		xi := data.Get(i)

		if xi < min {
			min = xi
			min_index = i
		}

		if xi > max {
			max = xi
			max_index = i
		}

		if math.IsNaN(xi) {
			min = xi
			max = xi
			min_index = i
			max_index = i
			break
		}
	}

	return
}
