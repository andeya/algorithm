// Copyright (c) Harri Rautila, 2012

// This file is part of github.com/hrautila/cvx package.
// It is free software, distributed under the terms of GNU Lesser General Public
// License Version 3, or any later version. See the COPYING tile included in this archive.

package sets

import (
	"fmt"
	"github.com/henrylee2cn/algorithm/matrix"
)

// FloatMatrixSet is a collection of named sets of float valued matrices.
type FloatMatrixSet struct {
	sets map[string][]*matrix.FloatMatrix
}

// Create new FloatMatrix collection with with empty named sets.
func FloatSetNew(names ...string) *FloatMatrixSet {
	sz := len(names)
	if sz == 0 {
		sz = 4
	}
	mcap := 2 * sz
	ms := new(FloatMatrixSet)
	ms.sets = make(map[string][]*matrix.FloatMatrix, mcap)
	for _, k := range names {
		ms.sets[k] = nil
	}
	return ms
}

// Create new FloatMatrix collection with with empty named sets.
func NewFloatSet(names ...string) *FloatMatrixSet {
	sz := len(names)
	if sz == 0 {
		sz = 4
	}
	mcap := 2 * sz
	ms := new(FloatMatrixSet)
	ms.sets = make(map[string][]*matrix.FloatMatrix, mcap)
	for _, k := range names {
		ms.sets[k] = nil
	}
	return ms
}

// Test if key is in the set.
func (M *FloatMatrixSet) Exists(name string) bool {
	_, ok := M.sets[name]
	return ok
}

// Get named set
func (M *FloatMatrixSet) At(name string) []*matrix.FloatMatrix {
	mset, _ := M.sets[name]
	return mset
}

// Set the contents of matrix set.
func (M *FloatMatrixSet) Set(key string, ms ...*matrix.FloatMatrix) {
	mset := make([]*matrix.FloatMatrix, 0)
	mset = append(mset, ms...)
	M.sets[key] = mset
}

// Append matrices to matrix set.
func (M *FloatMatrixSet) Append(key string, ms ...*matrix.FloatMatrix) {
	mset, ok := M.sets[key]
	if !ok {
		mset = make([]*matrix.FloatMatrix, 0)
	}
	mset = append(mset, ms...)
	M.sets[key] = mset
}

// Remove set.
func (M *FloatMatrixSet) Remove(key string) {
	delete(M.sets, key)
}

func (M *FloatMatrixSet) Keys() []string {
	s := make([]string, len(M.sets))
	for key := range M.sets {
		s = append(s, key)
	}
	return s
}

func (M *FloatMatrixSet) Copy() *FloatMatrixSet {
	Mnew := NewFloatSet(M.Keys()...)
	for _, key := range M.Keys() {
		mset := M.At(key)
		for _, m := range mset {
			Mnew.Append(key, m.Copy())
		}
	}
	return Mnew
}

func (M *FloatMatrixSet) Print() {
	for key := range M.sets {
		ms := M.sets[key]
		for i, m := range ms {
			fmt.Printf("** %s[%d] **\n%v\n", key, i, m)
		}
	}
}

func (M *FloatMatrixSet) PrintMatrix(format string) {
	for key := range M.sets {
		ms := M.sets[key]
		for i, m := range ms {
			fmt.Printf("** %s[%d] **\n%s\n", key, i, m.ToString(format))
		}
	}
}

// DimensionSet is a collection of named sets of sizes.
type DimensionSet struct {
	sets map[string][]int
}

// Create new dimension set with empty dimension info.
func DSetNew(names ...string) *DimensionSet {
	sz := len(names)
	if sz == 0 {
		sz = 4
	}
	mcap := 2 * sz

	ds := new(DimensionSet)
	ds.sets = make(map[string][]int, mcap)
	for _, k := range names {
		nset := make([]int, 0, 16)
		ds.sets[k] = nset
	}
	return ds
}

// Create new dimension set with empty dimension info.
func NewDimensionSet(names ...string) *DimensionSet {
	sz := len(names)
	if sz == 0 {
		sz = 4
	}
	mcap := 2 * sz

	ds := new(DimensionSet)
	ds.sets = make(map[string][]int, mcap)
	for _, k := range names {
		nset := make([]int, 0, 16)
		ds.sets[k] = nset
	}
	return ds
}

// Get named set
func (ds *DimensionSet) At(name string) []int {
	dset, _ := ds.sets[name]
	return dset
}

// Append sizes to dimension set key.
func (ds *DimensionSet) Append(key string, dims []int) {
	dset, ok := ds.sets[key]
	if !ok {
		dset = make([]int, 0, 2*len(dims))
	}
	for _, v := range dims {
		dset = append(dset, v)
	}
	ds.sets[key] = dset
}

// Append dimension key to dis.
func (ds *DimensionSet) Set(key string, dims []int) {
	dset := make([]int, 0, 2*len(dims))
	for _, v := range dims {
		dset = append(dset, v)
	}
	ds.sets[key] = dset
}

// Find maximum dimension in sets.
func (ds *DimensionSet) Max(keys ...string) int {
	mx := 0
	for _, key := range keys {
		dset := ds.sets[key]
		for _, n := range dset {
			if n > mx {
				mx = n
			}
		}
	}
	return mx
}

// Calculate sum over set of keys.
func (ds *DimensionSet) Sum(keys ...string) int {
	sz := 0
	for _, key := range keys {
		dset := ds.sets[key]
		for _, n := range dset {
			sz += n
		}
	}
	return sz
}

// Calculate sum of squared dimensions over set of keys.
func (ds *DimensionSet) SumSquared(keys ...string) int {
	sz := 0
	for _, key := range keys {
		dset := ds.sets[key]
		for _, n := range dset {
			sz += n * n
		}
	}
	return sz
}

// Calculate sum of packed dimensions over set of keys.
func (ds *DimensionSet) SumPacked(keys ...string) int {
	sz := 0
	for _, key := range keys {
		dset := ds.sets[key]
		for _, n := range dset {
			sz += n * (n + 1) / 2
		}
	}
	return sz
}

func (ds *DimensionSet) Printf() {
	for key := range ds.sets {
		fmt.Printf("set '%s': %v\n", key, ds.sets[key])
	}
}

// Local Variables:
// tab-width: 4
// End:
