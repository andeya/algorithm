// Copyright 2011,2012 Johann Höchtl. All rights reserved.
// Use of this source code is governed by a Modified BSD License
// that can be found in the LICENSE file.

// Automated tests for the cartconvert/osgb36 package
package osgb36

import (
	"fmt"
	"github.com/henrylee2cn/algorithm/cartconvert/cartconvert"
	"math"
	"testing"
)

// ## AOSGB36ToStruct
type oSGB36StringToStructParam struct {
	osgb36coord string
	prec        OSGB36prec
}

type oSGB36StringToStructTest struct {
	in  oSGB36StringToStructParam
	out *OSGB36Coord
}

var oSGB36StringToStructTestssuc = []oSGB36StringToStructTest{
	{
		oSGB36StringToStructParam{"NN000500", OSGB36Auto},
		NewOSGB36Coord("NN", 0, 500, 0, 3, OSGB36Auto),
	},
	{
		oSGB36StringToStructParam{"NN000510", OSGB36Auto},
		NewOSGB36Coord("NN", 0, 510, 0, 3, OSGB36Auto),
	},
	{
		oSGB36StringToStructParam{"NN1660071200", OSGB36Auto},
		NewOSGB36Coord("NN", 1660, 7120, 0, 4, OSGB36Auto),
	},
	{
		oSGB36StringToStructParam{"NN", OSGB36Auto},
		&OSGB36Coord{Zone: "NN", Easting: 0, Northing: 0, gridLen: 0},
	},
	{
		oSGB36StringToStructParam{"NN11", OSGB36Auto},
		&OSGB36Coord{Zone: "NN", Easting: 1, Northing: 1, gridLen: 1},
	},
	{
		oSGB36StringToStructParam{"NN1212", OSGB36Auto},
		&OSGB36Coord{Zone: "NN", Easting: 12, Northing: 12, gridLen: 2},
	},
	{
		oSGB36StringToStructParam{"NN123123", OSGB36Auto},
		&OSGB36Coord{Zone: "NN", Easting: 123, Northing: 123, gridLen: 3},
	},
	{
		oSGB36StringToStructParam{"NN12341234", OSGB36Auto},
		&OSGB36Coord{Zone: "NN", Easting: 1234, Northing: 1234, gridLen: 4},
	},
	{
		oSGB36StringToStructParam{"NN1234512345", OSGB36Auto},
		&OSGB36Coord{Zone: "NN", Easting: 12345, Northing: 12345, gridLen: 5},
	},
	{
		oSGB36StringToStructParam{"NN1234512345", OSGB36_2},
		&OSGB36Coord{Zone: "NN", Easting: 12, Northing: 12, gridLen: 2},
	},
	{
		oSGB36StringToStructParam{"NN1234512345", OSGB36_2},
		NewOSGB36Coord("NN", 12, 12, 0, 2, OSGB36Auto),
	},
	{
		oSGB36StringToStructParam{"NN166712", OSGB36_5},
		NewOSGB36Coord("NN", 1660, 7120, 0, 4, OSGB36_5),
	},
}

func osgb36equal(osgb1, osgb2 *OSGB36Coord) bool {
	p1 := fmt.Sprintf("%s", osgb1)
	p2 := fmt.Sprintf("%s", osgb2)
	return p1 == p2
}

func TestOSGB36StringToStruct(t *testing.T) {
	for cnt, test := range oSGB36StringToStructTestssuc {
		out, err := AOSGB36ToStruct(test.in.osgb36coord, test.in.prec)

		if err != nil {
			t.Errorf("AOSGB36ToStruct [%d]: Error: %s", cnt, err)
		} else {
			formattspec := "AOSGB36ToStruct [%d]: Expected %s, got %s"
			// fmt.Printf(formattspec+"\n", cnt, test.out, out)
			if !osgb36equal(test.out, out) {
				t.Errorf(formattspec, cnt, test.out, out)
			}
		}
	}
}

// ## OSGB36ToWGS84LatLong
type oSGB36ToWGS84LatLongTest struct {
	in  *OSGB36Coord
	out *cartconvert.PolarCoord
}

// Values taken from http://gridreferencefinder.com/
var oSGB36ToWGS84LatLongTests = []oSGB36ToWGS84LatLongTest{
	{
		&OSGB36Coord{Zone: "SE", Easting: 29793, Northing: 33798, gridLen: 5, El: cartconvert.Airy1830Ellipsoid},
		&cartconvert.PolarCoord{Latitude: 53.799638, Longitude: -1.5491515},
	},
	{
		&OSGB36Coord{Zone: "NN", Easting: 16600, Northing: 71200, gridLen: 3, El: cartconvert.Airy1830Ellipsoid},
		&cartconvert.PolarCoord{Latitude: 56.796557, Longitude: -5.0039304},
	},
	{
		&OSGB36Coord{Zone: "NN", Easting: 1, Northing: 7, gridLen: 5, El: cartconvert.Airy1830Ellipsoid},
		&cartconvert.PolarCoord{Latitude: 56.150684, Longitude: -5.2214376},
	},
	{
		&OSGB36Coord{Zone: "NN", Easting: 16600, Northing: 71200, gridLen: 5, El: cartconvert.Airy1830Ellipsoid},
		&cartconvert.PolarCoord{Latitude: 56.796088, Longitude: -5.0047120},
	},
	{
		&OSGB36Coord{Zone: "NN", Easting: 16650, Northing: 71250, gridLen: 5, El: cartconvert.Airy1830Ellipsoid},
		&cartconvert.PolarCoord{Latitude: 56.796557, Longitude: -5.0039304},
	},
	{
		&OSGB36Coord{Zone: "SV", Easting: 0, Northing: 0, gridLen: 0, El: cartconvert.Airy1830Ellipsoid},
		&cartconvert.PolarCoord{Latitude: 49.766809, Longitude: -7.5571598},
	},
}

func latlongequal(pcp1, pcp2 *cartconvert.PolarCoord) bool {
	// here we reduce precision to three digits. OSGB36 conversion using a helmert
	// transformation is only exacto to a maximum of 5m and for certain locations
	// may deviate from the true point to +/- 25m!
	pp1s := fmt.Sprintf("%.3g %.3g", pcp1.Latitude, pcp1.Longitude)
	pp2s := fmt.Sprintf("%.3g %.3g", pcp2.Latitude, pcp2.Longitude)

	return pp1s == pp2s
}

func TestOSGB36ToWGS84LatLong(t *testing.T) {
	for cnt, test := range oSGB36ToWGS84LatLongTests {

		out := OSGB36ToWGS84LatLong(test.in)

		if !latlongequal(test.out, out) {
			t.Errorf("OSGB36ToWGS84LatLong:%d [%s]: Expected %s, got %s", cnt, test.in, test.out, out)
		}
	}
}

// ## WGS84LatLongToOSGB36
type wGS84LatLongToOSGB36Test struct {
	in  *cartconvert.PolarCoord
	out *OSGB36Coord
}

// Values taken from http://gridreferencefinder.com/
var wGS84LatLongToOSGB36Tests = []wGS84LatLongToOSGB36Test{
	{
		&cartconvert.PolarCoord{Latitude: 53.799638, Longitude: -1.5491515},
		&OSGB36Coord{Zone: "SE", Easting: 29793, Northing: 33798, gridLen: 5, El: cartconvert.Airy1830Ellipsoid},
	},
	{
		&cartconvert.PolarCoord{Latitude: 56.796557, Longitude: -5.0039304},
		&OSGB36Coord{Zone: "NN", Easting: 16650, Northing: 71250, gridLen: 5, El: cartconvert.Airy1830Ellipsoid},
	},
}

func osgb36fuzzyequal(osgb1, osgb2 *OSGB36Coord) bool {
	if osgb1.Zone != osgb2.Zone {
		return false
	}
	// the helmert transformation between OSGB36 and WGS84 will not attain a higher accuracy but 5m
	// and exceed for certain locations +/- 25m. Therefore we accept a deviation of max +/- 15 as correct.
	return math.Sqrt(math.Pow(float64(int(osgb1.Easting-osgb2.Easting)), 2)+math.Pow(float64(int(osgb1.Northing-osgb2.Northing)), 2)) < 15.0
}

func TestWGS84LatLongToOSGB36(t *testing.T) {
	for cnt, test := range wGS84LatLongToOSGB36Tests {
		out, err := WGS84LatLongToOSGB36(test.in)
		if err != nil {
			t.Errorf("WGS84LatLongToOSGB36 [%d]: Error: %s", cnt, err)
		} else {
			if !osgb36fuzzyequal(test.out, out) {
				t.Errorf("WGS84LatLongToOSGB36:%d [%s]: Expected %s, got %s", cnt, test.in, test.out, out)
			}
		}
	}
}
