package units

import (
	"math"
	"time"
)

// float64 representation of distance in mm. Millimeters are used rather
// then meters to get an increased precision near 0.
type Distance float64

// FIXME: Rename to metre?
const (
	Nanometer           = 1e-6 * Millimeter
	Micrometer          = 1e-3 * Millimeter
	Millimeter Distance = 1
	Centimeter          = 1e1 * Millimeter
	Decimeter           = 1e2 * Millimeter
	Meter               = 1e3 * Millimeter
	Kilometer           = 1e6 * Millimeter
)

func (d Distance) Micrometers() float64 {
	return float64(d/Micrometer)
}

func (d Distance) Millimeters() float64 {
	return float64(d/Millimeter)
}

func (d Distance) Centimeters() float64 {
	return float64(d/Centimeter)
}

func (d Distance) Meters() float64 {
	return float64(d/Meter)
}

func (d Distance) Kilometers() float64 {
	return float64(d/Kilometer)
}

/* Operations */

func (d Distance) Abs() Distance {
	if(d < 0) {return -d}
	return d
}

func (d Distance) DivideWithDuration(t time.Duration) Velocity {
	return Velocity(float64(d)/float64(t)*float64(time.Second))
}

func Hypot(x, y Distance) Distance {
	return Distance(math.Hypot(float64(x), float64(y)))
}


// Backward compatibillity only!
//func DivideDistanceDuration(d Distance, t time.Duration) Velocity {
//	return d.DivideWithDuration(t)
//}

// Backward compatibillity only!
//func Abs(d Distance) Distance {
//	return d.Abs()
//}
