package stat

/* covariance.go
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

// takes a dataset and calculates the covariance
func covariance(data1, data2 Interface, mean1, mean2 float64) (res float64) {
	// calculate the sum of the squares
	n := data1.Len()
	for i := 0; i < n; i++ {
		delta1 := (data1.Get(i) - mean1)
		delta2 := (data2.Get(i) - mean2)
		res += (delta1*delta2 - res) / float64(i+1)
	}
	return
}

func CovarianceMean(data1, data2 Interface, mean1, mean2 float64) float64 {
	n := data1.Len()
	covariance := covariance(data1, data2, mean1, mean2)
	return covariance * float64(n) / float64(n-1)
}

func Covariance(data1, data2 Interface) float64 {
	mean1 := Mean(data1)
	mean2 := Mean(data2)

	return CovarianceMean(data1, data2, mean1, mean2)
}

/*
Correlation()
  Calculate Pearson correlation = cov(X, Y) / (sigma_X * sigma_Y)
This routine efficiently computes the correlation in one pass of the
data and makes use of the algorithm described in:

B. P. Welford, "Note on a Method for Calculating Corrected Sums of
Squares and Products", Technometrics, Vol 4, No 3, 1962.

This paper derives a numerically stable recurrence to compute a sum
of products

S = sum_{i=1..N} [ (x_i - mu_x) * (y_i - mu_y) ]

with the relation

S_n = S_{n-1} + ((n-1)/n) * (x_n - mu_x_{n-1}) * (y_n - mu_y_{n-1})
*/
func Correlation(data1, data2 Interface) (res float64) {
	var sum_xsq, sum_ysq, sum_cross float64
	n := data1.Len()

	//
	// Compute:
	// sum_xsq = Sum [ (x_i - mu_x)^2 ],
	// sum_ysq = Sum [ (y_i - mu_y)^2 ] and
	// sum_cross = Sum [ (x_i - mu_x) * (y_i - mu_y) ]
	// using the above relation from Welford's paper
	//

	mean_x := data1.Get(0)
	mean_y := data2.Get(0)

	for i := 1; i < n; i++ {
		ratio := float64(i) / float64(i+1)
		delta_x := data1.Get(i) - mean_x
		delta_y := data2.Get(i) - mean_y
		sum_xsq += delta_x * delta_x * ratio
		sum_ysq += delta_y * delta_y * ratio
		sum_cross += delta_x * delta_y * ratio
		mean_x += delta_x / float64(i+1)
		mean_y += delta_y / float64(i+1)
	}

	res = sum_cross / (math.Sqrt(sum_xsq) * math.Sqrt(sum_ysq))

	return
}
