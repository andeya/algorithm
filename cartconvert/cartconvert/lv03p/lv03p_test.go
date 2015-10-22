// Copyright 2011,2012 Johann Höchtl. All rights reserved.
// Use of this source code is governed by a Modified BSD License
// that can be found in the LICENSE file.

// Automated tests for the cartconvert/lv03p package
package lv03p

import (
	"fmt"
	"github.com/henrylee2cn/algorithm/cartconvert/cartconvert"
	"math"
	"testing"
)

type swissCoordRepresentation struct {
	in  SwissCoord
	out string
}

var swissCoordRepresentationTest = []swissCoordRepresentation{
	{
		SwissCoord{Easting: 235.5, Northing: 20.0, CoordType: LV03}, "y:235.5 x:20",
	},
}

func TestSwissCoordRepresentation(t *testing.T) {
	for cnt, test := range swissCoordRepresentationTest {
		out := test.in.String()

		if test.out != out {
			t.Errorf("TestSwissCoordRepresentation [%d]: Expected: %s, got: %s", cnt, test.out, out)
		}
	}
}

// ## ASwissCoordToStruct
type aSwissCoordToStructretparam struct {
	coord *SwissCoord
	err   error
}

func (val *aSwissCoordToStructretparam) String() (fs string) {

	if val.coord != nil {
		fs = val.coord.String()
	}

	if val.err != nil {
		fs += " " + val.err.Error()
	}
	return
}

type aSwissCoordToStruct struct {
	in  string
	out aSwissCoordToStructretparam
}

var aSwissCoordToStructTests = []aSwissCoordToStruct{
	{
		in: "x:25 y:34.3", out: aSwissCoordToStructretparam{coord: &SwissCoord{Easting: 34.3, Northing: 25, CoordType: LV03}, err: nil},
	},
	{
		in: "x:25.0 N:34.3", out: aSwissCoordToStructretparam{coord: nil, err: cartconvert.ErrSyntax},
	},
}

func aswisscoordtostructequal(coord1, coord2 aSwissCoordToStructretparam) bool {
	if coord1.coord != nil && coord2.coord != nil {
		return coord1.coord.String() == coord2.coord.String()
	}
	return coord1.err == coord2.err
}

func TestASwissCoordToStruct(t *testing.T) {
	for cnt, test := range aSwissCoordToStructTests {

		out, erro := ASwissCoordToStruct(test.in)
		retval := aSwissCoordToStructretparam{coord: out, err: erro}

		if !aswisscoordtostructequal(test.out, retval) {
			t.Errorf("ASwissCoordToStruct [%d]: expected %v, got %v", cnt, test.out, retval)
		}
	}
}

// ## SwissCoordToGRS80LatLong
type swissCoordToGRS80LatLongTest struct {
	in  *SwissCoord
	out *cartconvert.PolarCoord
}

// http://www.swisstopo.admin.ch/internet/swisstopo/de/home/apps/calc/navref.html
var swissCoordToGRS80LatLongTests = []swissCoordToGRS80LatLongTest{
	{
		&SwissCoord{Easting: 750536, Northing: 265013, CoordType: LV03, El: cartconvert.Bessel1841Ellipsoid},
		&cartconvert.PolarCoord{Latitude: 47.518605, Longitude: 9.437422},
	},
}

func latlongequal(pcp1, pcp2 *cartconvert.PolarCoord) bool {
	pp1s := fmt.Sprintf("%.4g %.4g", pcp1.Latitude, pcp1.Longitude)
	pp2s := fmt.Sprintf("%.4g %.4g", pcp2.Latitude, pcp2.Longitude)

	return pp1s == pp2s
}

func TestSwissCoordToGRS80LatLong(t *testing.T) {
	for cnt, test := range swissCoordToGRS80LatLongTests {

		out, _ := SwissCoordToGRS80LatLong(test.in)

		if !latlongequal(test.out, out) {
			t.Errorf("SwissCoordToGRS80LatLong [%d]: Expteced %s, got %s", cnt, test.out, out)
		}
	}
}

// ## GRS80LatLongToSwissCoord
type gRS80LatLongToSwissCoordParam struct {
	gc        *cartconvert.PolarCoord
	coordType SwissCoordType
}

type gRS80LatLongToSwissCoordTest struct {
	in  gRS80LatLongToSwissCoordParam
	out *SwissCoord
}

var gRS80LatLongToSwissCoordTests = []gRS80LatLongToSwissCoordTest{
	{
		gRS80LatLongToSwissCoordParam{
			gc:        &cartconvert.PolarCoord{Latitude: 47.518605, Longitude: 9.437422, El: cartconvert.GRS80Ellipsoid},
			coordType: LV03},
		NewSwissCoord(LV03, 750536, 265013, 0),
	},
}

func swisscoordfuzzyequal(c1, c2 *SwissCoord) bool {
	return math.Sqrt(math.Pow(c1.Easting-c2.Easting, 2)+math.Pow(c1.Northing-c2.Northing, 2)) < 2.0
}

func TestGRS80LatLongToSwissCoord(t *testing.T) {
	for cnt, test := range gRS80LatLongToSwissCoordTests {
		out, _ := GRS80LatLongToSwissCoord(test.in.gc, test.in.coordType)
		if swisscoordfuzzyequal(test.out, out) {
			t.Errorf("GRS80LatLongToSwissCoord [%d]: Expected %s, got %s", cnt, test.out, out)
		}
	}
}
