//
// GNU GSL Statistics library (v1.15, GPLv3) implemented in Go
//
package stat

/* types.go
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

//
// Interface is used throughout the package.
//
type Interface interface {
	Get(int) float64
	Len() int
}

//
// Float64Slice is a predifined float64 slice type
// which implements the Interface and Sort interfaces.
//
type Float64Slice []float64

func (f Float64Slice) Get(i int) float64 { return f[i] }
func (f Float64Slice) Len() int          { return len(f) }

// The following methods are *only* required for Sort and SortStrider

func (f Float64Slice) Less(i, j int) bool { return f[i] < f[j] }
func (f Float64Slice) Swap(i, j int)      { f[i], f[j] = f[j], f[i] }

//
// IntSlice is a predifined int64 slice type
// which implements the Interface and Sort interfaces.
//
type IntSlice []int64

func (f IntSlice) Get(i int) float64 { return float64(f[i]) }
func (f IntSlice) Len() int          { return len(f) }

// The following methods are *only* required for Sort and SortStrider

func (f IntSlice) Less(i, j int) bool { return f[i] < f[j] }
func (f IntSlice) Swap(i, j int)      { f[i], f[j] = f[j], f[i] }

//
// Strider strides over the data, for sampling purposes.
//
type Strider struct {
	Interface
	stride int
}

func NewStrider(data Interface, stride int) Strider {
	return Strider{data, stride}
}

func (p Strider) Get(i int) float64 {
	return p.Interface.Get(i * p.stride)
}

func (p Strider) Len() int {
	return p.Interface.Len() / p.stride
}

//
// Sort is Interface with sorting functionality.
//
type Sort interface {
	Interface
	Less(int, int) bool
	Swap(int, int)
}

//
// SortStrider strides over the data, for sampling purposes.
// It also has sorting functionality.
//
type SortStrider struct {
	Sort
	stride int
}

func NewSortStrider(data Sort, stride int) SortStrider {
	return SortStrider{data, stride}
}

func (p SortStrider) Get(i int) float64 {
	return p.Sort.Get(i * p.stride)
}

func (p SortStrider) Len() int {
	return p.Sort.Len() / p.stride
}

func (p SortStrider) Less(i, j int) bool {
	return p.Sort.Less(i*p.stride, j*p.stride)
}

func (p SortStrider) Swap(i, j int) {
	p.Sort.Swap(i*p.stride, j*p.stride)
}
