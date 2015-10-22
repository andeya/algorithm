package stat

/* wkurtosis.go
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

func WKurtosis(w, data Interface) float64 {
	wmean := WMean(w, data)
	wsd := WSdMean(w, data, wmean)
	return WKurtosisMeanSd(w, data, wmean, wsd)
}

// WKurtosisMean calculates the kurtosis of a dataset
func WKurtosisMeanSd(w, data Interface, wmean, wsd float64) float64 {

	var wavg, W float64
	Len := data.Len()

	for i := 0; i < Len; i++ {
		wi := w.Get(i)

		if wi > 0 {
			x := (data.Get(i) - wmean) / wsd
			W += wi
			wavg += (x*x*x*x - wavg) * (wi / W)
		}
	}

	return wavg - 3.0 // makes kurtosis zero for a Gaussian
}
