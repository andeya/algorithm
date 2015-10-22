package stat

/* 
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
	"testing"
)

var data = Float64Slice{1.0}
var data1 = Float64Slice{0.0, 1.0, 2.0, 3.0, 4.0, 5.0}
var data2 = make(Float64Slice, 10000)

var data_I = IntSlice{1}
var data1_I = IntSlice{0, 1, 2, 3, 4, 5}
var data2_I = make(IntSlice, 10000)

func BenchmarkMean(b *testing.B) {
	var f float64
	for i := 0; i < b.N; i++ {
		f = Mean(&data)
	}
	_ = f // make compiler happy
}

func BenchmarkMean1(b *testing.B) {
	var f float64
	for i := 0; i < b.N; i++ {
		f = Mean(&data1)
	}
	_ = f // make compiler happy
}

func BenchmarkMean2(b *testing.B) {
	var f float64
	for i := 0; i < b.N; i++ {
		f = Mean(&data2)
	}
	_ = f // make compiler happy
}

func BenchmarkMeanI(b *testing.B) {
	var f float64
	for i := 0; i < b.N; i++ {
		f = Mean(&data_I)
	}
	_ = f // make compiler happy
}

func BenchmarkMean1_I(b *testing.B) {
	var f float64
	for i := 0; i < b.N; i++ {
		f = Mean(&data1_I)
	}
	_ = f // make compiler happy
}

func BenchmarkMean2_I(b *testing.B) {
	var f float64
	for i := 0; i < b.N; i++ {
		f = Mean(&data2_I)
	}
	_ = f // make compiler happy
}
