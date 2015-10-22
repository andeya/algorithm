// Copyright 2011,2012 Johann Höchtl. All rights reserved.
// Use of this source code is governed by a Modified BSD License
// that can be found in the LICENSE file.

// This package provides functions to deal with conversion and transformations of coordinates
// of the OSGB36 datum, the UK National Grid.
//
// The conversion between WGS84 and OSGB36 uses a simple helmert transformation,
// which in the case of OSGB36 inconsistencies, may result in an accuracy not exceeding +/- 5m, and
// for cornercases worse than +/- 15m.
// If a higher accuracy is required, a set of helmert parameters must be used or the
// procedure described at http://www.ordnancesurvey.co.uk/gps/docs/Geomatics_world.pdf.
//
// For further info see http://gps.ordnancesurvey.co.uk/etrs89geo_natgrid.asp
package osgb36

import (
	"fmt"
	"github.com/henrylee2cn/algorithm/cartconvert/cartconvert"
	"math"
	"strconv"
	"strings"
)

// A OSGB36 coordinate is specified by zone, easting and northing.
type OSGB36Coord struct {
	Easting, Northing uint
	RelHeight         float64
	Zone              string
	El                *cartconvert.Ellipsoid
	gridLen           byte
}

// Controls formatting of an OSGB36 coordinate.
//
//    OSGB36Auto - for easting and northing the most compact representation will be found. Eg.: NN12006500 will be stored as NN1265
//    OSGB36Leave - the adverse of OSGB36Auto: Do not try to do any compacting on easting and northing
//    OSGB36_1 ... OSGB36_5 - set northing resp. easting to exactly 1 .. 5 digits. If the inputs are shorter than n,
//       they will be filled with '0'. N1256 with precision 5 will result in N1200056000. If the input is longer than n,
//       it will be truncated. N123567 with precision 2 will result in N1256. When shortening bearings, no rounding takes place.
//
// The precision of an OSGB36 coordinate can either be set explicitly; From meter - resolution (OSGB36_5)
// to the bare 100x100 km Zone OSGB36_1. OSGB36Min is a synonym for OSGB36_1, OSGB36_Max is a synonym for OSGB36_5.
type OSGB36prec byte

const (
	OSGB36_Min OSGB36prec = iota
	OSGB36_1
	OSGB36_2
	OSGB36_3
	OSGB36_4
	OSGB36_5
	OSGB36_Max  = OSGB36_5
	OSGB36Leave = 128
	OSGB36Auto  = OSGB36Leave + 1
)

// Canonical representation of a OSGB36 datum.
func (coord *OSGB36Coord) String() string {
	if coord.gridLen > 0 {
		return fmt.Sprintf("%s%0*d%0*d", coord.Zone, int(coord.gridLen), coord.Easting, int(coord.gridLen), coord.Northing)
	}
	return coord.Zone
}

// Parses a string representation of an OSGB36 coordinate datum into a OSGB36 coordinate struct. The literal
// can be specified as follows:
//    ZO EA NO
//    ZO EANO
//    ZOEANO
// with ZO the two letter zone specifier, EA easting and NO northing to the accuracy of a meter.
//
// For a description of prec see OSGB36prec.
//
// The reference ellipsoid of an OSGB36 will always be set to the Airy1830 ellipsoid.
//
// returns cartconvert.ErrSyntax if format is not understood
// returns cartconvert.ErrRange if values are outside the defined parameters for an OSGB36 bearing
func AOSGB36ToStruct(osgb36coord string, prec OSGB36prec) (*OSGB36Coord, error) {

	compact := strings.ToUpper(strings.TrimSpace(osgb36coord))
	var zone, enn string
	var east, north int
	var err error

L1:
	for _, item := range compact {
		switch {
		case item == ' ':
			continue L1
		case byte(item)-'0' >= 0 && byte(item)-'0' <= 9:
			enn += string(item)
		default:
			zone += string(item)
		}
	}

	zl := len(zone)
	if zl == 0 || zl > 2 {
		return nil, cartconvert.ErrSyntax
	}

	ennlen := byte(len(enn))
	if ennlen > 0 {
		if ennlen%2 > 0 {
			return nil, cartconvert.ErrRange
		}

		ennlen /= 2
		if ennlen > byte(OSGB36_Max) {
			return nil, cartconvert.ErrRange
		}

		east, err = strconv.Atoi(enn[:ennlen])
		if err != nil {
			return nil, err
		}
		north, err = strconv.Atoi(enn[ennlen:])
		if err != nil {
			return nil, err
		}
	}
	return NewOSGB36Coord(zone, uint(east), uint(north), 0, ennlen, prec), nil
}

// Returns northing and easting based on OSGB36 zone specifier relative to false northing and easting
func OSGB36ZoneToRefCoords(coord *OSGB36Coord) (easting, northing uint) {
	var l1, l2 uint

	// get numeric values of letter references, mapping A->0, B->1, C->2, etc:
	l1 = uint(coord.Zone[0] - 'A')
	if len(coord.Zone) > 1 {
		l2 = uint(coord.Zone[1] - 'A')
	}

	// shuffle down letters after 'I' since 'I' is not used in grid:
	if l1 > 7 {
		l1--
	}
	if l2 > 7 {
		l2--
	}

	// convert grid letters into 100km-square indexes from false origin (grid square SV):
	easting = (((l1-2)%5)*5 + l2%5) * 100000
	northing = ((19 - l1/5*5) - l2/5) * 100000

	// append numeric part of references to grid index:
	easting += coord.Easting
	northing += coord.Northing

	// if only the grid zone was specified, return the location at point.
	// for all other specified vales return location at middle of square
	fact := uint(math.Pow(10, float64(byte(OSGB36_Max)-coord.gridLen)))
	if fact < 10000 {
		easting += 5 * (fact / 10)
		northing += 5 * (fact / 10)
	}
	return
}

// Convert an OSGB36 coordinate value to a WGS84 based latitude and longitude coordinate.
//
// Important: A OSGB36 datum like NN11 will be internally expanded to NN1500015000 to point to the middle of the zone.
// For the point at NN1000010000 it is necessary to fully qualify northing and easting.
//
// Plain grid zone specifiers will NOT be shifted towards the middle of the square.
func OSGB36ToWGS84LatLong(coord *OSGB36Coord) *cartconvert.PolarCoord {

	easting, northing := OSGB36ZoneToRefCoords(coord)

	gc := cartconvert.InverseTransverseMercator(
		&cartconvert.GeoPoint{Y: float64(northing), X: float64(easting), El: coord.El},
		49,
		-2,
		0.9996012717,
		400000,
		-100000)

	cart := cartconvert.PolarToCartesian(gc)
	pt := cartconvert.HelmertWGS84ToOSGB36.InverseTransform(&cartconvert.Point3D{X: cart.X, Y: cart.Y, Z: cart.Z})

	return cartconvert.CartesianToPolar(&cartconvert.CartPoint{X: pt.X, Y: pt.Y, Z: pt.Z, El: cartconvert.WGS84Ellipsoid})
}

// Perform formating on an OSGB36 datum. For formatting see OSGB36prec.
func SanitizeOSGB36CoordToPrec(easting, northing *uint, inputprec byte, desiredprec OSGB36prec) byte {

	if *easting+*northing == 0 || inputprec == 0 {
		return byte(OSGB36_Min)
	}

	// inputpreclen := byte(max(int(uintlen(*northing)), int(uintlen(*easting))))

	if inputprec < byte(OSGB36_Max) {
		fact := uint(math.Pow(10, float64(byte(OSGB36_Max)-inputprec)))
		*easting *= fact
		*northing *= fact
	}

	switch desiredprec {
	case OSGB36Auto:
		northprec, eastprec := OSGB36_Max, OSGB36_Max
		easttmp, northtmp := *easting, *northing

		for easttmp%10 == 0 && eastprec > 0 {
			easttmp /= 10
			eastprec--
		}

		for northtmp%10 == 0 && northprec > 0 {
			northtmp /= 10
			northprec--
		}
		desiredprec = OSGB36prec(max(int(northprec), int(eastprec)))

		fallthrough

	case OSGB36_Min, OSGB36_1, OSGB36_2, OSGB36_3, OSGB36_4, OSGB36_5:
		if desiredprec < OSGB36_Max {
			fact := uint(math.Pow(10, float64(OSGB36_Max-desiredprec)))
			*easting /= fact
			*northing /= fact
		}
	case OSGB36Leave:
		desiredprec = OSGB36prec(inputprec)
	}
	return byte(desiredprec)
}

// Build OSGB36 coordinate from easting and northing relative to Grid. The parameter prec controls
// how the resulting OSGB36 coordinates are formated. See OSGB36prec.
//
// The function will return cartconvert.ErrRange if lat/long are not within the OSGB36 datum area.
func GridRefNumToLet(easting, northing uint, height float64, prec OSGB36prec) (*OSGB36Coord, error) {
	// get the 100km-grid indices
	easting100k := easting / 100000
	northing100k := northing / 100000

	if easting100k < 0 || easting100k > 6 || northing100k < 0 || northing100k > 12 {
		return nil, cartconvert.ErrRange
	}

	// translate those into numeric equivalents of the grid letters
	l1 := byte((19 - northing100k) - (19-northing100k)%5 + (easting100k+10)/5)
	l2 := byte((19-northing100k)*5%25 + easting100k%5)

	// compensate for skipped 'I' and build grid letter-pairs
	if l1 > 7 {
		l1++
	}
	if l2 > 7 {
		l2++
	}

	zone := string(l1+'A') + string(l2+'A')
	easting %= 100000
	northing %= 100000

	return NewOSGB36Coord(zone, easting, northing, 0, 5, OSGB36Leave), nil
}

// Transform a latitude / longitude coordinate datum into a OSGB36 coordinate.
//
// Important: The reference ellipsoid of the originating coordinate system will be assumed
// to be the WGS84Ellipsoid and will be set thereupon, regardless of the actually set reference ellipsoid.
func WGS84LatLongToOSGB36(gc *cartconvert.PolarCoord) (*OSGB36Coord, error) {
	// This sets the Ellipsoid to WGS84, regardless of the actual value set
	gc.El = cartconvert.WGS84Ellipsoid

	cart := cartconvert.PolarToCartesian(gc)
	pt := cartconvert.HelmertWGS84ToOSGB36.Transform(&cartconvert.Point3D{X: cart.X, Y: cart.Y, Z: cart.Z})
	polar := cartconvert.CartesianToPolar(&cartconvert.CartPoint{X: pt.X, Y: pt.Y, Z: pt.Z, El: cartconvert.Airy1830Ellipsoid})

	gp := cartconvert.DirectTransverseMercator(
		polar,
		49,
		-2,
		0.9996012717,
		400000,
		-100000)

	return GridRefNumToLet(uint(gp.X+0.5), uint(gp.Y+0.5), 0, OSGB36_Max)
}

func max(x, y int) int {
	if x > y {
		return x
	}
	return y
}

// Create a new OSGB36 coordinate from literals. The parameter prec plays an important role in how the literals
// are interpreted:
//    OSGB36Auto - for easting and northing the most compact representation will be found. Eg.: NN12006500 will be stored as NN1265
//    OSGB36Leave - the adverse of OSGB36Auto: Do not try to do any compacting on easting and northing
//    OSGB36_1 ... OSGB36_5 - set northing resp. easting to exactly 1 .. 5 digits. If the inputs are shorter than n,
//       they will be filled with '0'. N1256 with precision 5 will result in N1200056000. If the input is longer than n,
//       it will be truncated. N123567 with precision 2 will result in N1256. When shortening bearings, no rounding takes place.
//
func NewOSGB36Coord(Zone string, easting, northing uint, relheight float64, inputprec byte, desiredprec OSGB36prec) *OSGB36Coord {
	effbytes := SanitizeOSGB36CoordToPrec(&easting, &northing, inputprec, desiredprec)
	return &OSGB36Coord{Easting: easting, Northing: northing, RelHeight: relheight, Zone: Zone, gridLen: effbytes, El: cartconvert.Airy1830Ellipsoid}
}
