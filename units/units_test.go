package units

import (
	"testing"
	"time"
	"math"
)

// Used to check weather two floats are 'different'
const test_k = 0.0000001

func cmpf64(f1, f2 float64) bool {
	err := f2 -f1
	if err < 0 {
		err = -err
	}
	if err < test_k {
		return true
	}
	return false
}

func TestDistance(t *testing.T) {
	// Test some selected constants:
	if !cmpf64(float64(1e9 * Nanometer), float64(1 * Meter)) {
		t.Error("1e9 nm != 1 m")
	}
	if !cmpf64(float64(1000 * Millimeter), float64(1 * Meter)) {
		t.Error("1000 mm != 1 m")
	}
	if !cmpf64(float64(100 * Centimeter), float64(1 * Meter)) {
		t.Error("100 cm != 1 m")
	}
	if !cmpf64(float64(1000 * Meter), float64(1 * Kilometer)) {
		t.Error("10000 m != 1 km")
	}

	// Test f64 converter functions:
	var d Distance = 1 * Meter

	if !cmpf64(d.Kilometers(), 0.001) {
		t.Error("(1m).Kilometers() != 0.001")
	}
	if !cmpf64(d.Meters(), 1) {
		t.Error("(1m).Meters() != 1")
	}
	if !cmpf64(d.Centimeters(), 100) {
		t.Error("(1m).Centimeters() != 100")
	}
	if !cmpf64(d.Millimeters(), 1000) {
		t.Error("(1m).Millimeters() != 1000")
	}
	if !cmpf64(d.Micrometers(), 1e6) {
		t.Error("(1m).Micrometers() != 1e6")
	}
}

func TestVelocity(t *testing.T) {
	// TODO: Implement some tests
}

func TestAngles(t *testing.T) {
	a1 := 2 * math.Pi * Radian

	// Test all constants:
	a2 := 360 * Degree
	if !cmpf64(float64(a1), float64(a2)) {
		t.Error("2 Pi [rad] != 360 [deg]")
	}
	a2 = 400 * Gradian
	if !cmpf64(float64(a1), float64(a2)) {
		t.Error("2 Pi [rad] != 400 [grad]")
	}
	a2 = 2000 * math.Pi * Milliradian
	if !cmpf64(float64(a1), float64(a2)) {
		t.Error("2 Pi [rad] != 2000 Pi [mrad]")
	}

	// Test f64 functions:
	if !cmpf64(a1.Radians(), 2 * math.Pi) {
		t.Error("(2 Pi [rad]).Radians() != 2 Pi")
	}
	if !cmpf64(a1.Degrees(), 360) {
		t.Error("(2 Pi [rad]).Degrees() != 360")
	}
	if !cmpf64(a1.Gradians(), 400) {
		t.Error("(2 Pi [rad]).Gradians() != 400")
	}
}

func TestAngularVelocity(t *testing.T) {
	// TODO: Implement some tests
}

func TestUnitOperations(t *testing.T) {
	// velocity = distance / time
	v1 := Meter.DivideWithDuration(time.Second)
	v2 := MeterPerSecond
	if !cmpf64(float64(v1), float64(v2)) {
		t.Error("velocity != distance / time")
	}
	// distance = velocity * time
	d1 := MeterPerSecond.MultiplyWithDuration(time.Second)
	d2 := Meter
	if !cmpf64(float64(d1), float64(d2)) {
		t.Error("distance != velocity * time")
	}

	// angular velocity = angle / time
	r1 := Radian.DivideWithDuration(time.Second)
	r2 := RadianPerSecond
	if !cmpf64(float64(r1), float64(r2)) {
		t.Error("angular velocity != angle / time")
	}

	// angle = angular velocity * time
	a1 := RadianPerSecond.MultiplyWithDuration(time.Second)
	a2 := Radian
	if !cmpf64(float64(a1), float64(a2)) {
		t.Error("angle != angular velocity * time")
	}
}
