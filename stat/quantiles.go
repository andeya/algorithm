package stat

/* quantiles.go
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

// QuantileFromSortedData performs the quantile function, also called percent 
// point function or inverse cumulative distribution function, on the sorted data. 
// Note that the function doesn't check wheather the data is actually sorted.
func QuantileFromSortedData(sortedData Interface, f float64) (result float64) {
	n := sortedData.Len()
	index := f * float64(n-1)
	lhs := int(index)
	delta := index - float64(lhs)

	if n == 0 {
		return 0.0
	}

	if lhs == n-1 {
		result = sortedData.Get(lhs)
	} else {
		result = (1-delta)*sortedData.Get(lhs) + delta*sortedData.Get(lhs+1)
	}

	return
}
