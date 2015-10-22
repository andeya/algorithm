// Copyright 2011,2012 Johann Höchtl. All rights reserved.
// Use of this source code is governed by a Modified BSD License
// that can be found in the LICENSE file.

// This package provides a series of functions to deal with
// conversion and transformations of coordinates in the Swiss coordinate system.
//
// The Swiss coordinate system recently switched from lv03 to lv95, however the difference is
// generall below one meter and ignored by this package. For accuracy within 1cm the
//  FINELTRA-Transformation has to be applied.
//
// References:
//
// [DE]: http://www.swisstopo.admin.ch/internet/swisstopo/de/home/topics/survey/sys/refsys/switzerland.parsysrelated1.24280.downloadList.32633.DownloadFile.tmp/refsysd.pdf
package lv03p

import (
	"fmt"
	"github.com/henrylee2cn/algorithm/cartconvert/cartconvert"
	"strconv"
	"strings"
)

// Coordinate type of Switzerland. Only affects string representation but not accuracy (The two systems diverge by about 1m)
type SwissCoordType byte

const (
	LV03 SwissCoordType = iota
	LV95
)

// A coordinate in Switzerland is specified by easting (right-value, x), and northing (height-value, y)
type SwissCoord struct {
	Easting, Northing, RelHeight float64
	CoordType                    SwissCoordType
	El                           *cartconvert.Ellipsoid
}

var coordliterals = [][]string{{"y:", " x:"}, {"E:", " N:"}}

// Canonical representation of a SwissCoord-value
func (bc *SwissCoord) String() (fs string) {

	var next float64

	if bc == nil {
		return
	}
	for i := 0; i < 2; i++ {
		fs += coordliterals[bc.CoordType][i]
		switch i {
		case 0:
			next = bc.Easting
		case 1:
			next = bc.Northing
		}

		tmp := fmt.Sprintf("%f", next)
		n := len(tmp)
		for n > 0 && tmp[n-1] == '0' {
			n--
		}
		if n > 0 && tmp[n-1] == '.' {
			n--
		}
		fs = fs + tmp[:n]
	}
	return
}

// Parses a string representation of a LV++ coordinate into a struct holding a SwissCoord coordinate value.
// The reference ellipsoid of Swisscoord datum is always the GRS80 ellipsoid.
func ASwissCoordToStruct(coord string) (*SwissCoord, error) {

	compact := strings.ToUpper(strings.TrimSpace(coord))
	var rights, heights string
	var coordType, oldcoordType SwissCoordType
	var right, height float64
	var err error

L1:
	for i, index := 0, 0; i < 2; i++ {
		index = strings.Index(compact, " ")
		if index == -1 {
			index = len(compact)
		}

		switch compact[:2] {
		case "X:":
			coordType = LV03
			heights = compact[2:index]
		case "Y:":
			coordType = LV03
			rights = compact[2:index]
		case "E:":
			coordType = LV95
			rights = compact[2:index]
		case "N:":
			coordType = LV95
			heights = compact[2:index]
		default:
			err = cartconvert.ErrSyntax
			break L1
		}

		if oldcoordType != coordType {
			err = cartconvert.ErrSyntax
			break L1
		}

		if i == 1 {
			break L1
		}
		compact = compact[index+len(" "):]
		compact = strings.TrimLeft(compact, " ")
		oldcoordType = coordType
	}

	if err == nil {

		right, err = strconv.ParseFloat(rights, 64)
		if err == nil {

			height, err = strconv.ParseFloat(heights, 64)
			if err == nil {
				return &SwissCoord{Easting: right, Northing: height, CoordType: coordType, El: cartconvert.Bessel1841Ellipsoid}, nil
			}
		}
	}

	return nil, err
}

// Transform a Swiss coordinate value to a GRS80 based latitude and longitude coordinate. Function returns
// cartconvert.ErrRange, if the swiss coordinate type is not one of LV03 or LV95
func SwissCoordToGRS80LatLong(coord *SwissCoord) (*cartconvert.PolarCoord, error) {

	var fn, fe float64

	switch coord.CoordType {
	case LV03:
		fe = 600000
		fn = 200000
	case LV95:
		fe = -2600000
		fn = -1200000
	default:
		return nil, cartconvert.ErrRange
	}

	gc := cartconvert.InverseTransverseMercator(
		&cartconvert.GeoPoint{Y: coord.Northing, X: coord.Easting, El: coord.El},
		46.952406, // lat0
		7.439583,  // long0
		1,
		fe, // fe
		fn) // fn

	cart := cartconvert.PolarToCartesian(gc)
	// According to literature, the Granit87 parameters shall not be used in favour of
	// higher accuracy of the following shift values

	// pt := cartconvert.HelmertLV03ToWGS84Granit87.Transform(&cartconvert.Point3D{X: cart.X, Y: cart.Y, Z: cart.Z})
	pt := &cartconvert.Point3D{X: cart.X + 674.374, Y: cart.Y + 15.056, Z: cart.Z + 405.346}

	return cartconvert.CartesianToPolar(&cartconvert.CartPoint{X: pt.X, Y: pt.Y, Z: pt.Z, El: cartconvert.GRS80Ellipsoid}), nil
}

// Transform a latitude / longitude coordinate datum into a Swiss coordinate. Function returns
// cartconvert.ErrRange, if the coordinate type is not set.
//
// Important: The reference ellipsoid of the originating coordinate system will be assumed
// to be the GRS80Ellipsoid and will be set thereupon, regardless of the actually set reference ellipsoid.
func GRS80LatLongToSwissCoord(gc *cartconvert.PolarCoord, coordType SwissCoordType) (*SwissCoord, error) {

	var fn, fe float64

	// This sets the Ellipsoid to GRS80, regardless of the actual value set
	gc.El = cartconvert.GRS80Ellipsoid

	cart := cartconvert.PolarToCartesian(gc)
	// According to literature, the Granit87 parameters shall not be used in favour of
	// higher accuracy of the following shift values

	// pt := cartconvert.HelmertWGS84ToMGI.Transform(&cartconvert.Point3D{X: cart.X, Y: cart.Y, Z: cart.Z})
	pt := &cartconvert.Point3D{X: cart.X - 674.374, Y: cart.Y - 15.056, Z: cart.Z - 405.346}
	polar := cartconvert.CartesianToPolar(&cartconvert.CartPoint{X: pt.X, Y: pt.Y, Z: pt.Z, El: cartconvert.Bessel1841Ellipsoid})

	switch coordType {
	case LV03:
		fe = 600000
		fn = 200000
	case LV95:
		fe = -2600000
		fn = -1200000
	default:
		return nil, cartconvert.ErrRange
	}

	gp := cartconvert.DirectTransverseMercator(
		polar,
		46.952406, // lat0
		7.439583,  // long0
		1,
		fe, // fe
		fn) // fn

	return &SwissCoord{CoordType: coordType, Northing: gp.Y, Easting: gp.X, El: gp.El}, nil
}

func NewSwissCoord(CoordType SwissCoordType, Easting, Northing, RelHeight float64) *SwissCoord {
	return &SwissCoord{Easting: Easting, Northing: Northing, RelHeight: RelHeight, CoordType: CoordType, El: cartconvert.Bessel1841Ellipsoid}
}
