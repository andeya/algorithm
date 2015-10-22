// Copyright 2011,2012 Johann Höchtl. All rights reserved.
// Use of this source code is governed by a Modified BSD License
// that can be found in the LICENSE file.

// Automated tests for the cartconvert package
package cartconvert

import (
	"fmt"
	"math"
	"testing"
)

// ## TestDirectTransverseMercator
type directtransversemercatorParam struct {
	pc                         *PolarCoord
	lat0, long0, scale, fe, fn float64
}

type directtransversemercatorTest struct {
	in  directtransversemercatorParam
	out *GeoPoint
}

var directtransversemercatorTests = []directtransversemercatorTest{
	{
		directtransversemercatorParam{
			&PolarCoord{Latitude: 50.5, Longitude: 0.5, El: Airy1830Ellipsoid},
			49,
			-2,
			0.9996012717,
			400000,
			-100000,
		},
		&GeoPoint{X: 577274.9838, Y: 69740.4923},
	},
}

func geopointequal(gp1, gp2 *GeoPoint) bool {
	gp1s := fmt.Sprintf("%.4f %.4f", gp1.X, gp1.Y)
	gp2s := fmt.Sprintf("%.4f %.4f", gp2.X, gp2.Y)

	return gp1s == gp2s
}

func TestDirectTransverseMercator(t *testing.T) {
	for cnt, test := range directtransversemercatorTests {
		out := DirectTransverseMercator(test.in.pc, test.in.lat0, test.in.long0, test.in.scale, test.in.fe, test.in.fn)
		if !geopointequal(test.out, out) {
			t.Errorf("DirectTransverseMercator [%d]: Expected %v, got %v", cnt, test.out, out)
		}
	}
}

// ## InverseTransverseMercator
type inversetransversemercatorParam struct {
	pt                         *GeoPoint
	lat0, long0, scale, fe, fn float64
}

type inversetransversemercatorTest struct {
	in  inversetransversemercatorParam
	out *PolarCoord
}

var inversetransversemercatorTests = []inversetransversemercatorTest{
	{
		inversetransversemercatorParam{
			&GeoPoint{X: 577274.99, Y: 69740.5, El: Airy1830Ellipsoid},
			49,
			-2,
			0.9996012717,
			400000,
			-100000,
		},
		&PolarCoord{Latitude: 50.5, Longitude: 0.5},
	},
}

func latlongequal(pcp1, pcp2 *PolarCoord) bool {
	pp1s := fmt.Sprintf("%.5f %.5f", pcp1.Latitude, pcp1.Longitude)
	pp2s := fmt.Sprintf("%.5f %.5f", pcp2.Latitude, pcp2.Longitude)

	return pp1s == pp2s
}

func TestInverseTransverseMercator(t *testing.T) {
	for _, test := range inversetransversemercatorTests {
		out := InverseTransverseMercator(test.in.pt, test.in.lat0, test.in.long0, test.in.scale, test.in.fe, test.in.fn)
		if !latlongequal(test.out, out) {
			t.Errorf("InverseTransverseMercator")
		}
	}
}

// ## ADegMMSSToNum
type degMMSSToNumTest struct {
	in  string
	out float64
}

var degMMSSToNumTests = []degMMSSToNumTest{
	{" N 30 °56 ' 34.45 ''", 30.942903},
	{" N30 ° 56 ' 34.45'' ", 30.942903},
	{" N 170°56'34.45''", 170.942903},
	{"  170 °", 170.0},
	{" - 359 ° 30'", -359.500000},
	{" - 180°30'30''", -180.508333},
	{" - 180°0'30''", -180.008333},
	{" - 180°00'30''", -180.008333},
	{" - 180°0'0.5''", -180.000139},
	{" - 180°30'", -180.5},
}

func floatequal(f1, f2 float64) bool {
	sf1 := fmt.Sprintf("%f", f1)
	sf2 := fmt.Sprintf("%f", f2)

	return sf1 == sf2
}

func TestDegMMSSToNum(t *testing.T) {
	for index, test := range degMMSSToNumTests {
		out, err := ADegMMSSToNum(test.in)

		if err != nil {
			t.Errorf("Error@ADegMMSSToNum [%d]: %s", index, err)
		}

		if !floatequal(test.out, out) {
			t.Errorf("ADegMMSSToNum [%d]: expected %f, got %f", index, test.out, out)
		}
	}
}

// ## ADegCommaToNum
type degCommaToNumTest struct {
	in  string
	out float64
}

var degCommaToNumTests = []degMMSSToNumTest{
	{" - 179.50°", -179.5},
	{" - 179  ° ", -179.0},
	{" S 50.50  °", -50.5},
}

func TestDegCommaToNum(t *testing.T) {
	for index, test := range degCommaToNumTests {
		out, err := ADegCommaToNum(test.in)

		if err != nil {
			t.Error(err)
		}

		if !floatequal(test.out, out) {
			t.Errorf("ADegCommaToNum %d: expected %f, got %f", index, test.out, out)
		}
	}
}

// ## LatLongToUTM
type aLatLongToUTMTest struct {
	in  *PolarCoord
	out *UTMCoord
}

var aLatLongToUTMTests = []aLatLongToUTMTest{
	{
		in:  &PolarCoord{Latitude: -33.922667, Longitude: 18.416689},
		out: &UTMCoord{Zone: "34H", Easting: 261190.0, Northing: 6243413.0},
	},
	{
		in:  &PolarCoord{Latitude: 43.642567, Longitude: -79.387139},
		out: &UTMCoord{Zone: "17T", Easting: 630084.0, Northing: 4833439.0},
	},
}

func utmabrequal(utm1, utm2 *UTMCoord) bool {
	u1 := fmt.Sprintf("%s %.0f %.0f", utm1.Zone, utm1.Easting, utm1.Northing)
	u2 := fmt.Sprintf("%s %.0f %.0f", utm2.Zone, utm2.Easting, utm2.Northing)

	return u1 == u2
}

func TestLatLongToUTM(t *testing.T) {
	for _, test := range aLatLongToUTMTests {
		out := LatLongToUTM(test.in)
		if !utmabrequal(test.out, out) {
			t.Errorf("LatLongToUTM")
		}
	}
}

// ## AUTMToStruct
type aUTMToStructTestParam struct {
	utmcoord string
	el       *Ellipsoid
}

type aUTMToStructTest struct {
	in  aUTMToStructTestParam
	out *UTMCoord
}

var aUTMToStructTests = []aUTMToStructTest{
	{
		aUTMToStructTestParam{" 17T   630084.31    4833438.548 ", WGS84Ellipsoid},
		&UTMCoord{Zone: "17T", Easting: 630084.31, Northing: 4833438.548},
	},
	{
		aUTMToStructTestParam{" 34H   261190    6243413 ", WGS84Ellipsoid},
		&UTMCoord{Zone: "34H", Easting: 261190.0, Northing: 6243413.0},
	},
}

func utmequal(utm1, utm2 *UTMCoord) bool {
	u1 := fmt.Sprintf("%s %f %f", utm1.Zone, utm1.Easting, utm1.Northing)
	u2 := fmt.Sprintf("%s %f %f", utm2.Zone, utm2.Easting, utm2.Northing)

	return u1 == u2
}

func TestAUTMToStruct(t *testing.T) {
	for _, test := range aUTMToStructTests {
		out, err := AUTMToStruct(test.in.utmcoord, test.in.el)

		if err != nil {
			t.Error(err)
		}

		if !utmequal(test.out, out) {
			t.Error("UTMToStruct")
		}
	}
}

// ## UTMToLatLong
type uTMToLatLongTest struct {
	in  *UTMCoord
	out *PolarCoord
}

var uTMToLatLongTests = []uTMToLatLongTest{
	{
		&UTMCoord{Zone: "17T", Easting: 630084.31, Northing: 4833438.548},
		&PolarCoord{Latitude: 43.642567, Longitude: -79.387139},
	},
	{
		&UTMCoord{Zone: "34H", Easting: 261190.0, Northing: 6243413.0},
		&PolarCoord{Latitude: -33.922665, Longitude: 18.416688},
	},
}

func polarequal(pc1, pc2 *PolarCoord) bool {
	p1 := fmt.Sprintf("%f %f", pc1.Latitude, pc1.Longitude)
	p2 := fmt.Sprintf("%f %f", pc2.Latitude, pc2.Longitude)

	return p1 == p2
}

func TestUTMToLatLong(t *testing.T) {
	for _, test := range uTMToLatLongTests {
		out, err := UTMToLatLong(test.in)

		if err != nil {
			t.Error(err)
		}

		if !polarequal(test.out, out) {
			t.Error("UTMToLatLong")
		}
	}
}

// ## ADegMMSSToPolar
type aDegMMSSToPolarParam struct {
	Northing, Easting string
	Height            float64
	El                *Ellipsoid
}

type aDegMMSSToPolarTest struct {
	in  aDegMMSSToPolarParam
	out *PolarCoord
}

var aDegMMSSToPolarTests = []aDegMMSSToPolarTest{
	{
		aDegMMSSToPolarParam{"S33° 55' 21.6''", "E18° 25' 0.08''", 0, WGS84Ellipsoid},
		&PolarCoord{Latitude: -33.922667, Longitude: 18.416689},
	},
	{
		aDegMMSSToPolarParam{"N43°38'33.24''", "W79°23'13.7''", 0, WGS84Ellipsoid},
		&PolarCoord{Latitude: 43.642567, Longitude: -79.387139},
	},
}

func TestADegMMSSToPolar(t *testing.T) {
	for _, test := range aDegMMSSToPolarTests {
		out, err := ADegMMSSToPolar(test.in.Northing, test.in.Easting, test.in.Height, test.in.El)

		if err != nil {
			t.Error(err)
		}

		if !latlongequal(test.out, out) {
			t.Error("ADegMMSSToPolar")
		}
	}
}

// ## PolarToCartesian
type polarToCartesianTest struct {
	in  *PolarCoord
	out *CartPoint
}

var polarToCartesianTests = []polarToCartesianTest{
	{
		&PolarCoord{Latitude: 47.567, Longitude: 14.243, El: WGS84Ellipsoid},
		&CartPoint{Y: 1060748.224300, X: 4178845.984047, Z: 4684527.101880},
	},
}

func cartequal(cp1, cp2 *CartPoint) bool {
	p1 := fmt.Sprintf("%f %f %f", cp1.Y, cp1.X, cp1.Z)
	p2 := fmt.Sprintf("%f %f %f", cp2.Y, cp2.X, cp2.Z)

	return p1 == p2
}

func TestPolarToCartesian(t *testing.T) {
	for index, test := range polarToCartesianTests {
		out := PolarToCartesian(test.in)
		if !cartequal(test.out, out) {
			t.Errorf("PolarToCartesian [%d]: expected (%f %f %f), got (%f %f %f)",
				index, test.out.Y, test.out.X, test.out.Z, out.Y, out.X, out.Z)
		}
	}
}

// ## CartesianToPolar
type cartesianToPolarTest struct {
	in  *CartPoint
	out *PolarCoord
}

var cartesianToPolarTests = []cartesianToPolarTest{
	{
		&CartPoint{Y: 1060748.224300, X: 4178845.984047, Z: 4684527.101880, El: WGS84Ellipsoid},
		&PolarCoord{Latitude: 47.567, Longitude: 14.243},
	},
}

func TestCartesianToPolar(t *testing.T) {
	for index, test := range cartesianToPolarTests {
		out := CartesianToPolar(test.in)
		if !latlongequal(test.out, out) {
			t.Errorf("CartesianToPolar [%d]: expected (%f %f), got (%f %f)",
				index, test.out.Latitude, test.out.Longitude, out.Latitude, out.Longitude)
		}
	}
}

// ## Helmert
type helmertTest struct {
	in  *Point3D
	out *Point3D
}

var helmertTests = []helmertTest{
	{
		&Point3D{Y: 4178845.984047, X: 1060748.2243, Z: 4684527.10188},
		&Point3D{Y: 4178655.486121, X: 1060094.493596, Z: 4684148.315584},
	},
}

func point3dequal(p1, p2 *Point3D) bool {
	pt1 := fmt.Sprintf("%f %f %f", p1.X, p1.Y, p1.Z)
	pt2 := fmt.Sprintf("%f %f %f", p2.X, p2.Y, p2.Z)

	return pt1 == pt2
}

// the inversion of a helmert is not totaly accurate; let's assume an accuracy of 2 meters
// (may even be more)
func helmertfuzzyequal(p1, p2 *Point3D) bool {
	return math.Sqrt(math.Pow(p1.X-p2.X, 2)+math.Pow(p1.Y-p2.Y, 2)+math.Pow(p1.Z-p2.Z, 2)) < 2.0
}

func TestHelmertTransform(t *testing.T) {
	for index, test := range helmertTests {
		out := HelmertWGS84ToMGI.Transform(test.in)
		if !point3dequal(test.out, out) {
			t.Errorf("HelmertTransform [%d]: expected (%f %f %f), got (%f %f %f)",
				index, test.out.X, test.out.Y, test.out.Z, out.X, out.Y, out.Z)
		}
	}
}

func TestHelmertInverseTransform(t *testing.T) {
	for index, test := range helmertTests {
		out := HelmertWGS84ToMGI.InverseTransform(test.out)
		if !helmertfuzzyequal(test.in, out) {
			t.Errorf("HelmertInverseTransform [%d]: expected (%f %f %f), got (%f %f %f)",
				index, test.in.X, test.in.Y, test.in.Z, out.X, out.Y, out.Z)
		}
	}
}

// ## GeoHashToLatLong
type geoHashToLatLongTest struct {
	in  string
	out *PolarCoord
}

var geoHashToLatLongTests = []geoHashToLatLongTest{
	{"ezs42",
		&PolarCoord{Latitude: 42.6, Longitude: -5.6},
	},
	{"u4pruydqqvj",
		&PolarCoord{Latitude: 57.64911, Longitude: 10.40744},
	},
	{"vn0p",
		&PolarCoord{Latitude: 80.0, Longitude: 45.0},
	},
	{"60p0",
		&PolarCoord{Latitude: -45.0, Longitude: -80.0},
	},
}

func TestGeoHashToLatLong(t *testing.T) {
	for index, test := range geoHashToLatLongTests {
		out, err := GeoHashToLatLong(test.in, nil)

		if err != nil {
			t.Error(err)
		}

		if !latlongequal(test.out, out) {
			t.Errorf("GeoHashToLatLong [%d]: expected (%f, %f), got (%f, %f)", index, test.out.Latitude, test.out.Longitude, out.Latitude, out.Longitude)
		}
	}
}

// LatLongToGeoHash - uses the inverted input of TestGeoHashToLatLong
func TestLatLongToGeoHash(t *testing.T) {
	for index, test := range geoHashToLatLongTests {
		out := LatLongToGeoHash(test.out)

		if test.in != out {
			t.Errorf("LatLongToGeoHash [%d]: expected %s, got %s", index, test.in, out)
		}
	}
}

// ## LatLongToString
type latlongToStringTest struct {
	in        *PolarCoord
	lat, long string
}

var latlongToStringTests = []latlongToStringTest{
	{&PolarCoord{Latitude: 140.0, Longitude: 40.0}, "N 140°", "E 40°"},
	{&PolarCoord{Latitude: -140.0, Longitude: -40.0}, "S 140°", "W 40°"},
	{&PolarCoord{Latitude: 140.5, Longitude: -40.5}, "N 140°30'", "W 40°30'"},
	{&PolarCoord{Latitude: 140.005, Longitude: -40.505}, "N 140°0'18''", "W 40°30'18''"},
}

func TestLatLongToString(t *testing.T) {
	for index, test := range latlongToStringTests {
		lat, long := LatLongToString(test.in, LLFdms)

		if !(test.lat == lat && test.long == long) {
			t.Errorf("LatLongToString [%d]: expected %s, %s, got %s, %s", index, test.lat, test.long, lat, long)
		}
	}
}
