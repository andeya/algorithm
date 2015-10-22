package stat_test

/* example_test.go
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
	"fmt"
	"github.com/henrylee2cn/algorithm/stat"
	"sort"
)

func ExampleMean() {
	data := stat.Float64Slice{17.2, 18.1, 16.5, 18.3, 12.6}
	mean := stat.Mean(data)
	fmt.Printf("The sample mean is %g", mean)
	// Output: The sample mean is 16.54
}

func ExampleVariance() {
	data := stat.Float64Slice{17.2, 18.1, 16.5, 18.3, 12.6}
	variance := stat.Variance(data)
	fmt.Printf("The estimated variance is %.4f", variance)
	// Output: The estimated variance is 5.3730
}

func ExampleMax() {
	data := stat.Float64Slice{17.2, 18.1, 16.5, 18.3, 12.6}
	largest, index := stat.Max(data)
	fmt.Printf("The largest value is %g and the index is %d", largest, index)
	// Output: The largest value is 18.3 and the index is 3
}

func ExampleMin() {
	data := stat.Float64Slice{17.2, 18.1, 16.5, 18.3, 12.6}
	smallest, index := stat.Min(data)
	fmt.Printf("The smallest value is %g and the index is %d", smallest, index)
	// Output: The smallest value is 12.6 and the index is 4
}

func ExampleMedianFromSortedData() {
	data := stat.Float64Slice{17.2, 18.1, 16.5, 18.3, 12.6}
	sort.Sort(data)
	median := stat.MedianFromSortedData(data)
	fmt.Printf("Sorted dataset: %v\n", data)
	fmt.Printf("The median is %g\n", median)
	// Output:
	// Sorted dataset: [12.6 16.5 17.2 18.1 18.3]
	// The median is 17.2
}

func ExampleQuantileFromSortedData() {
	data := stat.Float64Slice{17.2, 18.1, 16.5, 18.3, 12.6}
	sort.Sort(data)
	upperq := stat.QuantileFromSortedData(data, 0.75)
	lowerq := stat.QuantileFromSortedData(data, 0.25)

	fmt.Printf("Sorted dataset: %v\n", data)
	fmt.Printf("The upper quartile is %g\n", upperq)
	fmt.Printf("The lower quartile is %g\n", lowerq)
	// Output:
	// Sorted dataset: [12.6 16.5 17.2 18.1 18.3]
	// The upper quartile is 18.1
	// The lower quartile is 16.5
}

func ExampleStrider() {
	data := stat.Float64Slice{
		.0421, .0941, .1064, .0242,
		.1331, .0773, .0243, .0815,
		.1186, .0356, .0728, .0999,
		.0614, .0479}
	strider := stat.NewStrider(data, 4)
	for i := 0; i < strider.Len(); i++ {
		fmt.Println(strider.Get(i))
	}
	// Output:
	// 0.0421
	// 0.1331
	// 0.1186
}

func ExampleNewStrider() {
	data := stat.Float64Slice{
		.0421, .0941, .1064, .0242,
		.1331, .0773, .0243, .0815,
		.1186, .0356, .0728, .0999,
		.0614, .0479}
	strider := stat.NewStrider(data, 4)
	fmt.Printf("mean data is %.4f\n", stat.Mean(data))
	fmt.Printf("mean strider is %.4f\n", stat.Mean(strider))
	// Output:
	// mean data is 0.0728
	// mean strider is 0.0979
}
