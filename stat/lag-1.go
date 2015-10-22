package stat

/* lag-1.go
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

func Lag1Autocorrelation(data Interface) float64 {
	mean := Mean(data)
	return Lag1AutocorrelationMean(data, mean)
}

func Lag1AutocorrelationMean(data Interface, mean float64) float64 {
	var r1, q float64
	Len := data.Len()
	v := (data.Get(0) - mean) * (data.Get(0) - mean)

	for i := 1; i < Len; i++ {
		delta0 := data.Get(i-1) - mean
		delta1 := data.Get(i) - mean
		q += (delta0*delta1 - q) / float64(i+1)
		v += (delta1*delta1 - v) / float64(i+1)
	}

	r1 = q / v

	return r1
}
