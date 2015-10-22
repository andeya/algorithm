package stat

/* variance.go
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

func _variance(data Interface, mean float64) (variance float64) {
	Len := data.Len()

	// calculate the sum of the squares
	for i := 0; i < Len; i++ {
		delta := data.Get(i) - mean
		// TODO: long double for variance... How to implement in Go?
		variance += ((delta * delta) - variance) / float64(i+1)
	}

	return
}

func VarianceWithFixedMean(data Interface, mean float64) float64 {
	return _variance(data, mean)
}

func SdWithFixedMean(data Interface, mean float64) float64 {
	variance := _variance(data, mean)
	return math.Sqrt(variance)
}

func VarianceMean(data Interface, mean float64) float64 {
	variance := _variance(data, mean)
	return variance * float64(data.Len()) / float64(data.Len()-1)
}

func SdMean(data Interface, mean float64) float64 {
	variance := _variance(data, mean)
	return math.Sqrt(variance * float64(data.Len()) / float64(data.Len()-1))
}

func Variance(data Interface) float64 {
	mean := Mean(data)
	return VarianceMean(data, mean)
}

func Sd(data Interface) float64 {
	mean := Mean(data)
	return SdMean(data, mean)
}

// TssMean takes a dataset and finds the sum of squares about the mean
func TssMean(data Interface, mean float64) (res float64) {

	// find the sum of the squares
	n := data.Len()
	for i := 0; i < n; i++ {
		delta := data.Get(i) - mean
		res += delta * delta
	}

	return
}

func Tss(data Interface) float64 {
	mean := Mean(data)
	return TssMean(data, mean)
}
