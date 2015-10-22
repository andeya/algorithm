// Copyright 2011,2012 Johann Höchtl. All rights reserved.
// Use of this source code is governed by a Modified BSD License
// that can be found in the LICENSE file.

// Automated tests for the cartconvert/bmn (Bundesmeldenetz) package
package bmn

import (
	"fmt"
	"github.com/henrylee2cn/algorithm/cartconvert/cartconvert"
	"testing"
)

// ## BMNStringToStruct
type bMNStringToStructTest struct {
	in  string
	out *BMNCoord
}

var bMNStringToStructTests = []bMNStringToStructTest{
	{
		"M31 592269 272290",
		&BMNCoord{Meridian: BMNM31, Right: 592269.0, Height: 272290.0},
	},
}

func bmnequal(bmn1, bmn2 *BMNCoord) bool {
	p1 := fmt.Sprintf("%s", bmn1)
	p2 := fmt.Sprintf("%s", bmn2)

	return p1 == p2
}

func TestBMNStringToStruct(t *testing.T) {
	for _, test := range bMNStringToStructTests {
		out, err := ABMNToStruct(test.in)

		if err != nil {
			t.Error(err)
		}

		if !bmnequal(test.out, out) {
			t.Error("BMNStringToStruct")
		}
	}
}

// ## BMNToWGS84LatLong
type bMNToWGS84LatLongTest struct {
	in  *BMNCoord
	out *cartconvert.PolarCoord
}

func bMNStringToStructHelper(coord string) (bmncoord *BMNCoord) {
	bmncoord, _ = ABMNToStruct(coord)
	return
}

var bMNToWGS84LatLongTests = []bMNToWGS84LatLongTest{
	{
		NewBMNCoord(BMNM28, 592270.0, 272290, 0),
		&cartconvert.PolarCoord{Latitude: 47.439212, Longitude: 16.197434},
	},
	{ // TODO: Ist das möglich??
		bMNStringToStructHelper("M34 592269 272290"),
		&cartconvert.PolarCoord{Latitude: 47.570299, Longitude: 14.236188},
	},
	{
		bMNStringToStructHelper("M34 703168 374510"),
		&cartconvert.PolarCoord{Latitude: 48.507001, Longitude: 15.698748},
	},
}

func latlongequal(pcp1, pcp2 *cartconvert.PolarCoord) bool {
	pp1s := fmt.Sprintf("%.5g %.5g", pcp1.Latitude, pcp1.Longitude)
	pp2s := fmt.Sprintf("%.5g %.5g", pcp2.Latitude, pcp2.Longitude)

	return pp1s == pp2s
}

func TestBMNToWGS84LatLong(t *testing.T) {
	for _, test := range bMNToWGS84LatLongTests {

		out, _ := BMNToWGS84LatLong(test.in)

		if !latlongequal(test.out, out) {
			t.Error("BMNToWGS84LatLong")
		}
	}
}

// ## WGS84LatLongToBMN
type wGS84LatLongToBMNParam struct {
	gc       *cartconvert.PolarCoord
	meridian BMNMeridian
}

type wGS84LatLongToBMNTest struct {
	in  wGS84LatLongToBMNParam
	out *BMNCoord
}

var wGS84LatLongToBMNTests = []wGS84LatLongToBMNTest{
	{
		wGS84LatLongToBMNParam{
			gc:       &cartconvert.PolarCoord{Latitude: 47.570299, Longitude: 14.236188, El: cartconvert.WGS84Ellipsoid},
			meridian: BMNM34},
		NewBMNCoord(BMNM34, 592269, 272290.05, 0),
	},
	{
		wGS84LatLongToBMNParam{
			gc:       &cartconvert.PolarCoord{Latitude: 48.507001, Longitude: 15.698748, El: cartconvert.WGS84Ellipsoid},
			meridian: BMNZoneDet},
		NewBMNCoord(BMNM34, 703168, 374510, 0),
	},
}

func TestWGS84LatLongToBMN(t *testing.T) {
	for index, test := range wGS84LatLongToBMNTests {
		out, _ := WGS84LatLongToBMN(test.in.gc, test.in.meridian)
		if !bmnequal(test.out, out) {
			t.Errorf("WGS84LatLongToBMN [%d]: expected %s, got %s", index, test.out, out)
		}
	}
}
