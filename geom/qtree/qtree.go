// Copyright 2012 The geom Authors. All rights reserved.
// Use of this source code is governed by a license that
// can be found in the LICENSE file.

package qtree

import (
	"fmt"
	"github.com/henrylee2cn/algorithm/geom"
)

var Debug = false
var Indent = 0

func dbg(format string, args ...interface{}) {
	if Debug {
		buf := make([]byte, Indent)
		for i := range buf {
			buf[i] = ' '
		}
		fmt.Printf(string(buf)+format+"\n", args...)
	}
}

type Config struct {
	SplitCount     int     "Exclusive upper bound before splitting a node"
	SplitSizeRatio float64 "Elements must be small enough to split"
	Height         int     "Max depth of search tree"
}

func ConfigDefault() (cfg Config) {
	cfg.SplitCount = 5
	cfg.SplitSizeRatio = 0.5
	cfg.Height = 10
	return
}

type Item interface {
	geom.Bounded
	Equals(oi interface{}) bool
}

type Tree struct {
	height int
	cfg    Config

	Count int

	//a promise that nothing added will be outside this bounds
	UpperBounds geom.Rect

	//smallest rect that contains all current items
	Bounds geom.Rect

	Partition geom.Coord

	Subtrees [4]*Tree

	Elements    map[Item]bool
	BigElements map[Item]bool
}

func New(cfg Config, bounds geom.Rect) (me *Tree) {
	me = &Tree{
		cfg:         cfg,
		UpperBounds: bounds,
		Partition:   bounds.Min.Plus(bounds.Max).Times(0.5),
	}
	me.Bounds = geom.NilRect()
	return
}

func (me *Tree) IsBig(bounds geom.Rect) bool {
	return bounds.Width() >= me.cfg.SplitSizeRatio*me.UpperBounds.Width() ||
		bounds.Height() >= me.cfg.SplitSizeRatio*me.UpperBounds.Height()
}

func (me *Tree) insertSubTrees(element Item) (inserted bool) {
	for i, t := range me.Subtrees {
		if t == nil {
			subbounds := me.UpperBounds
			switch i {
			case 0:
				subbounds.Min.X = me.Partition.X
				subbounds.Min.Y = me.Partition.Y
			case 1:
				subbounds.Min.X = me.Partition.X
				subbounds.Max.Y = me.Partition.Y
			case 2:
				subbounds.Max.X = me.Partition.X
				subbounds.Min.Y = me.Partition.Y
			case 3:
				subbounds.Max.X = me.Partition.X
				subbounds.Max.Y = me.Partition.Y
			}
			cfg := me.cfg
			cfg.Height--
			t = New(cfg, subbounds)
			me.Subtrees[i] = t
		}

		inserted = inserted || t.Insert(element)
	}

	return
}

func (me *Tree) Find(element Item) (found Item, ok bool) {
	if !geom.RectsIntersect(element.Bounds(), me.UpperBounds) {
		return
	}
	if me.IsBig(element.Bounds()) && me.BigElements != nil {
		for elem := range me.BigElements {
			if element.Equals(elem) {
				found = elem
				ok = true
				return
			}
		}
		return
	}
	if me.Elements != nil {
		for elem := range me.Elements {
			if element.Equals(elem) {
				found = elem
				ok = true
				return
			}
		}
		return
	}
	for _, t := range me.Subtrees {
		if t == nil {
			continue
		}
		found, ok = t.Find(element)
		if ok {
			return
		}
	}
	return
}

func (me *Tree) FindOrInsert(element Item) (found Item, inserted bool) {
	defer func(bounds geom.Rect) {
		if inserted {
			me.Bounds.ExpandToContainRect(bounds)
			me.Count++
		}
	}(element.Bounds())

	if !geom.RectsIntersect(element.Bounds(), me.UpperBounds) {
		dbg("doesn't belong in %v", me.UpperBounds)
		return
	}

	if me.IsBig(element.Bounds()) {
		if me.BigElements == nil {
			me.BigElements = make(map[Item]bool)
		}
		for elem := range me.BigElements {
			//if geom.RectsEqual(elem.Bounds(), element.Bounds()) {
			if element.Equals(elem) {
				found = elem
				inserted = false
				return
			}
		}
		me.BigElements[element] = true
		found = element
		inserted = true
		return
	}

	if me.cfg.Height == 0 {
		if me.Elements == nil {
			me.Elements = make(map[Item]bool)
		}
	} else {
		if me.Elements != nil && len(me.Elements) == me.cfg.SplitCount {
			for elem := range me.Elements {
				me.insertSubTrees(elem)
			}
			me.Elements = nil
		}
	}

	if me.Subtrees[0] != nil {
		dbg("looking through subtrees")
		for _, t := range me.Subtrees {
			if t == nil {
				continue
			}
			foundInSubtree, insertedInSubtree := t.FindOrInsert(element)
			if foundInSubtree != nil {
				// if found in the subtree, all subtrees should agree here
				found = foundInSubtree
				inserted = insertedInSubtree
			}
		}

		return
	}

	if me.Elements != nil {
		dbg("looking through Element")
		for elem := range me.Elements {
			//if geom.RectsEqual(elem.Bounds(), element.Bounds()) {
			if element.Equals(elem) {
				found = elem
				inserted = false
				return
			}
		}
	} else {
		me.Elements = make(map[Item]bool)
	}

	me.Elements[element] = true
	found = element
	inserted = true
	return
}

func (me *Tree) Insert(element Item) (inserted bool) {
	str := ""
	if geom.RectsIntersect(element.Bounds(), me.UpperBounds) {
		str = "*"
	}
	dbg("inserting in %v%s", me.Bounds, str)

	if !geom.RectsIntersect(me.UpperBounds, element.Bounds()) {
		return
	}

	defer func(bounds geom.Rect) {
		me.Bounds.ExpandToContainRect(bounds)
		me.Count++
	}(element.Bounds())

	inserted = true

	//if this element is too big, stop here
	if me.IsBig(element.Bounds()) {
		if me.BigElements == nil {
			me.BigElements = make(map[Item]bool)
		}
		me.BigElements[element] = true
		return
	}

	//if we're at the bottom, stop here
	if me.cfg.Height == 0 {
		if me.Elements == nil {
			me.Elements = make(map[Item]bool)
		}
		me.Elements[element] = true
		return
	}

	//if we've got enough at this level, break into subtrees
	if me.Elements != nil && len(me.Elements) == me.cfg.SplitCount {
		for elem := range me.Elements {
			me.insertSubTrees(elem)
		}
		me.Elements = nil
	}

	//if we already have subtrees, insert into them
	if me.Subtrees[0] != nil {
		me.insertSubTrees(element)
		return
	}

	//no subtrees, stop here
	if me.Elements == nil {
		me.Elements = make(map[Item]bool)
	}
	me.Elements[element] = true

	return
}

func (me *Tree) Remove(element Item) (removed bool) {
	dbg("removing %v", element)
	if Debug {
		println(element)
	}
	Indent++
	defer func() { Indent-- }()

	if !geom.RectsIntersect(me.UpperBounds, element.Bounds()) {
		dbg("out of bounds")
		return
	}

	defer func() {
		if removed {
			me.Count = 0
			me.Bounds = geom.NilRect()
			if me.BigElements != nil {
				me.Count += len(me.BigElements)
				for elem := range me.BigElements {
					me.Bounds.ExpandToContainRect(elem.Bounds())
				}
			}
			if me.Elements != nil {
				me.Count += len(me.Elements)
				for elem := range me.Elements {
					me.Bounds.ExpandToContainRect(elem.Bounds())
				}
			}
			if me.Subtrees[0] != nil {
				for _, t := range me.Subtrees {
					if t.Count != 0 {
						me.Count += t.Count
						me.Bounds.ExpandToContainRect(t.Bounds)
					}
				}
			}
		}
	}()

	dbg("BigElements: %v", me.BigElements)
	dbg("Elements: %v", me.Elements)

	if me.BigElements != nil {
		for elem := range me.BigElements {
			if element.Equals(elem) {
				delete(me.BigElements, elem)
				removed = true
			}
		}
	}
	if me.Elements != nil {
		for elem := range me.Elements {
			if element.Equals(elem) {
				delete(me.Elements, elem)
				removed = true
			}
		}
	}
	for _, t := range me.Subtrees {
		if t == nil {
			continue
		}
		if t.Remove(element) {
			removed = true
		}
	}

	return
}

func (me *Tree) Iterate() (rch <-chan Item) {
	ch := make(chan Item, me.Count)
	visited := make(map[Item]bool)
	me.iterate(ch, visited)
	close(ch)
	return ch
}
func (me *Tree) iterate(sch chan<- Item, visited map[Item]bool) {
	if me.BigElements != nil {
		for e := range me.BigElements {
			if !visited[e] {
				sch <- e
				visited[e] = true
			}
		}
	}
	if me.Elements != nil {
		for e := range me.Elements {
			if !visited[e] {
				sch <- e
				visited[e] = true
			}
		}
	}
	for _, t := range me.Subtrees {
		if t != nil {
			t.iterate(sch, visited)
		}
	}
}

func (me *Tree) Items() (rch <-chan Item) {
	col := make(map[Item]bool, me.Count)
	me.Enumerate(col)
	ch := make(chan Item, len(col))
	for i := range col {
		ch <- i
	}
	close(ch)
	rch = ch
	return
}

func (me *Tree) Enumerate(collection map[Item]bool) {
	if me.Elements != nil {
		for elem := range me.Elements {
			collection[elem] = true
		}
	}

	if me.BigElements != nil {
		for elem := range me.BigElements {
			collection[elem] = true
		}
	}

	for _, t := range me.Subtrees {
		if t != nil {
			t.Enumerate(collection)
		}
	}
}

func (me *Tree) Do(foo func(x Item)) {
	if me.Elements != nil {
		for elem := range me.Elements {
			foo(elem)
		}
	}

	if me.BigElements != nil {
		for elem := range me.BigElements {
			foo(elem)
		}
	}

	for _, t := range me.Subtrees {
		if t != nil {
			t.Do(foo)
		}
	}
}

func (me *Tree) CollectInside(bounds geom.Rect, collection map[Item]bool) {
	if !geom.RectsIntersect(bounds, me.UpperBounds) {
		return
	}
	if me.BigElements != nil {
		for elem := range me.BigElements {
			if bounds.ContainsRect(elem.Bounds()) {
				collection[elem] = true
			}
		}
	}
	if me.Elements != nil {
		for elem := range me.Elements {
			if bounds.ContainsRect(elem.Bounds()) {
				collection[elem] = true
			}
		}
	}

	for _, t := range me.Subtrees {
		if t == nil {
			continue
		}
		t.CollectInside(bounds, collection)
	}

	return
}

func (me *Tree) RemoveInside(bounds geom.Rect, collection map[Item]bool) {
	if !geom.RectsIntersect(bounds, me.UpperBounds) {
		return
	}
	if me.BigElements != nil {
		for elem := range me.BigElements {
			if bounds.ContainsRect(elem.Bounds()) {
				delete(me.BigElements, elem)
				if collection != nil {
					collection[elem] = true
				}
			}
		}
	}
	if me.Elements != nil {
		for elem := range me.Elements {
			if bounds.ContainsRect(elem.Bounds()) {
				delete(me.Elements, elem)
				if collection != nil {
					collection[elem] = true
				}
			}
		}
	}

	for _, t := range me.Subtrees {
		if t == nil {
			continue
		}
		t.RemoveInside(bounds, collection)
	}

	return
}

func (me *Tree) CollectIntersect(bounds geom.Rect, collection map[Item]bool) {
	if !geom.RectsIntersect(bounds, me.UpperBounds) {
		return
	}
	if me.BigElements != nil {
		for elem := range me.BigElements {
			if geom.RectsIntersect(elem.Bounds(), bounds) {
				collection[elem] = true
			}
		}
	}

	if me.Elements != nil {
		for elem := range me.Elements {
			if geom.RectsIntersect(bounds, elem.Bounds()) {
				collection[elem] = true
			}
		}
	}

	for _, t := range me.Subtrees {
		if t == nil {
			continue
		}
		t.CollectIntersect(bounds, collection)
	}

	return
}

func (me *Tree) RemoveIntersect(bounds geom.Rect, collection map[Item]bool) {
	if !geom.RectsIntersect(bounds, me.UpperBounds) {
		return
	}
	if me.BigElements != nil {
		for elem := range me.BigElements {
			if geom.RectsIntersect(elem.Bounds(), bounds) {
				delete(me.BigElements, elem)
				if collection != nil {
					collection[elem] = true
				}
			}
		}
	}

	if me.Elements != nil {
		for elem := range me.Elements {
			if geom.RectsIntersect(bounds, elem.Bounds()) {
				delete(me.Elements, elem)
				if collection != nil {
					collection[elem] = true
				}
			}
		}
	}

	for _, t := range me.Subtrees {
		if t == nil {
			continue
		}
		t.RemoveIntersect(bounds, collection)
	}

	return
}

func (me *Tree) Size() int {
	return me.Count
}

func (me *Tree) String() string {
	str := "[]"
	if me.Subtrees[0] != nil {
		str = fmt.Sprintf("%v", me.Subtrees)
	}
	return fmt.Sprintf("QTree{%v, %v, %v, %s}", me.UpperBounds, me.Elements, me.BigElements, str)
}
