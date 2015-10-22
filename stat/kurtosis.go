package stat

/* kurtosis.go
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

func Kurtosis(data Interface) float64 {
	mean := Mean(data)
	est_sd := SdMean(data, mean)
	return KurtosisMainSd(data, mean, est_sd)
}

func KurtosisMainSd(data Interface, mean, sd float64) float64 {
	var avg, kurtosis float64
	Len := data.Len()

	// calculate the fourth moment the deviations, normalized by the sd

	// we use a recurrence relation to stable update a running value so
	// there aren't any large sums that can overflow

	for i := 0; i < Len; i++ {
		x := (data.Get(i) - mean) / sd
		avg += (x*x*x*x - avg) / float64(i+1)
	}

	kurtosis = avg - 3.0 // makes kurtosis zero for a Gaussian

	return kurtosis
}
