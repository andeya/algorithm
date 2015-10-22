// Author: slowpoke <proxypoke at lavabit dot com>
// Repository: https://gist.github.com/proxypoke/vector
//
// This program is free software under the non-terms
// of the Anti-License. Do whatever the fuck you want.

package vector

import (
	"strconv"
)

// For missmatched dimensions.
type DimError struct {
	Dim_a, Dim_b uint
}

type (
	// For out-of-range indexes.
	IndexError uint
	// For cross products where either ndim != 3.
	CrossError DimError
)

func (e DimError) Error() string {
	return "Missmatched dimensions: " +
		strconv.Itoa(int(e.Dim_b)) +
		" != " +
		strconv.Itoa(int(e.Dim_b))
}

func (e IndexError) Error() string {
	return "Index out of range: " + strconv.Itoa(int(e))
}

func (e CrossError) Error() string {
	return "Invalid dimensions: " +
		strconv.Itoa(int(e.Dim_a)) +
		", " +
		strconv.Itoa(int(e.Dim_b)) +
		" (must be 3)"
}
