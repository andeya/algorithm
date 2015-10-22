package stat

/* skew.go
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

func Skew(data Interface) float64 {
	mean := Mean(data)
	sd := SdMean(data, mean)
	return SkewMeanSd(data, mean, sd)
}

// SkewMeanSd calculates the skewness of a dataset
func SkewMeanSd(data Interface, mean, sd float64) (skew float64) {
	Len := data.Len()

	for i := 0; i < Len; i++ {
		x := (data.Get(i) - mean) / sd
		skew += (x*x*x - skew) / float64(i+1)
	}
	return
}
