// Copyright 2011,2012 Johann Höchtl. All rights reserved.
// Use of this source code is governed by a Modified BSD License
// that can be found in the LICENSE file.

// This package provides a series of functions to deal with
// conversion and transformations of coordinates in the Datum Austria
//
// and here specifically of the Bundesmeldenetz, the former federal cartographic datum
// of Austria. The Bundesmeldenetz is already widely replaced by UTM coordinates but much legacy
// data is still encoded in BMN coordinates. Unlike UTM, the BMN uses the Bessel reference ellipsoid
// and uses lat0 at Hierro (canary islands), which makes transformations tedious.
// For more information see
//
// [DE]: http://www.topsoft.at/pstrainer/entwicklung/algorithm/karto/oek/austria_oek.htm#bmn
// [EN]: http://www.asprs.org/resources/grids/03-2004-austria.pdf
package bmn

import (
	"fmt"
	"github.com/henrylee2cn/algorithm/cartconvert/cartconvert"
	"strconv"
	"strings"
)

// Meridian Coordinates of the Bundesmeldenetz, three values describing false easting and false northing.
// The meridian specification of BMN plays the same role as the zone specifier of UTM.
type BMNMeridian byte

const (
	BMNZoneDet BMNMeridian = iota
	BMNM28
	BMNM31
	BMNM34
)

func (bm BMNMeridian) String() (rep string) {
	switch bm {
	case BMNM28:
		rep = "M28"
	case BMNM31:
		rep = "M31"
	case BMNM34:
		rep = "M34"
	case BMNZoneDet:
		rep = "autodetect"
	default:
		rep = "#unknown"
	}
	return
}

// A BMN coordinate is specified by right-value (easting), height-value (northing)
// and the meridian stripe, 28°, 31° or 34° West of Hierro
type BMNCoord struct {
	Right, Height, RelHeight float64
	Meridian                 BMNMeridian
	El                       *cartconvert.Ellipsoid
}

// Canonical representation of a BMN-value
func (bc *BMNCoord) String() (fs string) {

	fs = bc.Meridian.String()
	var next float64

	for i := 0; i < 2; i++ {
		fs += " "
		switch i {
		case 0:
			next = bc.Right
		case 1:
			next = bc.Height
		}

		tmp := fmt.Sprintf("%.0f", next)
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

// Parses a string representation of a BMN-Coordinate into a struct holding a BMN coordinate value.
// The reference ellipsoid of BMN coordinates is always the Bessel ellipsoid.
func ABMNToStruct(bmncoord string) (*BMNCoord, error) {

	compact := strings.ToUpper(strings.TrimSpace(bmncoord))
	var rights, heights string
	var meridian BMNMeridian
	var right, height float64
	var err error

L1:
	for i, index := 0, 0; i < 3; i++ {
		index = strings.Index(compact, " ")
		if index == -1 {
			index = len(compact)
		}
		switch i {
		case 0:
			switch compact[:index] {
			case "M28":
				meridian = BMNM28
			case "M31":
				meridian = BMNM31
			case "M34":
				meridian = BMNM34
			default:
				err = cartconvert.ErrSyntax
				break L1
			}
		case 1:
			rights = compact[:index]
		case 2:
			heights = compact[:index]
			break L1
		}
		compact = compact[index+len(" "):]
		compact = strings.TrimLeft(compact, " ")
	}

	if err == nil {

		right, err = strconv.ParseFloat(rights, 64)
		if err == nil {

			height, err = strconv.ParseFloat(heights, 64)
			if err == nil {

				return &BMNCoord{Right: right, Height: height, Meridian: meridian, El: cartconvert.Bessel1841MGIEllipsoid}, nil
			}
		}
	}

	return nil, err
}

// Transform a BMN coordinate value to a WGS84 based latitude and longitude coordinate. Function returns
// cartconvert.ErrRange, if the meridian stripe of the bmn-coordinate is not set
func BMNToWGS84LatLong(bmncoord *BMNCoord) (*cartconvert.PolarCoord, error) {

	var long0, fe float64

	switch bmncoord.Meridian {
	case BMNM28:
		long0 = 10.0 + 20.0/60.0
		fe = 150000
	case BMNM31:
		long0 = 13.0 + 20.0/60.0
		fe = 450000
	case BMNM34:
		long0 = 16.0 + 20.0/60.0
		fe = 750000
	default:
		return nil, cartconvert.ErrRange
	}

	gc := cartconvert.InverseTransverseMercator(
		&cartconvert.GeoPoint{Y: bmncoord.Height, X: bmncoord.Right, El: bmncoord.El},
		0,
		long0,
		1,
		fe,
		-5000000)

	cart := cartconvert.PolarToCartesian(gc)
	pt := cartconvert.HelmertWGS84ToMGI.InverseTransform(&cartconvert.Point3D{X: cart.X, Y: cart.Y, Z: cart.Z})

	return cartconvert.CartesianToPolar(&cartconvert.CartPoint{X: pt.X, Y: pt.Y, Z: pt.Z, El: cartconvert.WGS84Ellipsoid}), nil
}

// Transform a latitude / longitude coordinate datum into a BMN coordinate. Function returns
// cartconvert.ErrRange, if the meridian stripe of the bmn-coordinate is not set.
//
// Important: The reference ellipsoid of the originating coordinate system will be assumed
// to be the WGS84Ellipsoid and will be set thereupon, regardless of the actually set reference ellipsoid.
func WGS84LatLongToBMN(gc *cartconvert.PolarCoord, meridian BMNMeridian) (*BMNCoord, error) {

	var long0, fe float64

	// This sets the Ellipsoid to WGS84, regardless of the actual value set
	gc.El = cartconvert.WGS84Ellipsoid

	cart := cartconvert.PolarToCartesian(gc)
	pt := cartconvert.HelmertWGS84ToMGI.Transform(&cartconvert.Point3D{X: cart.X, Y: cart.Y, Z: cart.Z})
	polar := cartconvert.CartesianToPolar(&cartconvert.CartPoint{X: pt.X, Y: pt.Y, Z: pt.Z, El: cartconvert.Bessel1841MGIEllipsoid})

	// Determine meridian stripe based on longitude
	if meridian == BMNZoneDet {
		switch {
		case 11.0+0.5/6*10 >= polar.Longitude && polar.Longitude >= 8.0+0.5/6*10:
			meridian = BMNM28
		case 14.0+0.5/6*10 >= polar.Longitude && polar.Longitude >= 11.0+0.5/6*10:
			meridian = BMNM31
		case 17.0+0.5/6*10 >= polar.Longitude && polar.Longitude >= 14.0+0.5/6*10:
			meridian = BMNM34
		}
	}

	switch meridian {
	case BMNM28:
		long0 = 10.0 + 20.0/60.0
		fe = 150000
	case BMNM31:
		long0 = 13.0 + 20.0/60.0
		fe = 450000
	case BMNM34:
		long0 = 16.0 + 20.0/60.0
		fe = 750000
	default:
		return nil, cartconvert.ErrRange
	}

	gp := cartconvert.DirectTransverseMercator(
		polar,
		0,
		long0,
		1,
		fe,
		-5000000)

	return &BMNCoord{Meridian: meridian, Height: gp.Y, Right: gp.X, El: gp.El}, nil
}

func NewBMNCoord(Meridian BMNMeridian, Right, Height, RelHeight float64) *BMNCoord {
	return &BMNCoord{Right: Right, Height: Height, RelHeight: RelHeight, Meridian: Meridian, El: cartconvert.Bessel1841MGIEllipsoid}
}
