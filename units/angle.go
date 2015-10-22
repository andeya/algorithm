package units

import (
	"math"
	"time"
)

// An Angle represented in radians. An integer representation would give
// better predictability (under flow), but almost all angles needs the
// math-library and atan2...
type Angle float64

const (
	Milliradian  Angle = 1
	Radian       Angle = 1e3
	Degree             = 2*math.Pi*Radian/360.0
	Gradian            = 2*math.Pi*Radian/400.0
)

func (a Angle) Radians() float64 {
	return float64(a/Radian)
}

func (a Angle) Degrees() float64 {
	return float64(a/Degree)
}

func (a Angle) Gradians() float64 {
	return float64(a/Gradian)
}

func (a Angle) DivideWithDuration(t time.Duration) AngularVelocity {
	return AngularVelocity(float64(a)/float64(t)*float64(time.Second))
}

func (a Angle) Abs() Angle {
	if(a < 0) {return -a}
	return a
}

func Atan2(y, x Distance) Angle {
	return Angle(math.Atan2(float64(y), float64(x)))*Radian
}

// Normalize the angle to be between -pi and +pi degrees.
func (a Angle) Normalize() Angle {
	return Angle(NormAngle(a.Radians()))*Radian
}

// Normalize an angle to range +/- math.Pi. Assumes angle is in radians
func NormAngle(angle float64) float64 {
	angle = math.Mod(angle, 2.0*math.Pi)
	if angle > math.Pi {
		angle -= 2*math.Pi
	} else if angle < -math.Pi {
		angle += 2*math.Pi
	}
	return angle
}

// Backwards compatibility only!
//func DivideAngleDuration(a Angle, t time.Duration) AngularVelocity {
//	return AngularVelocity(float64(a)/float64(t)*float64(time.Second))
//}

// Backwards compatibility only!
//func (a Angle) NormAngle() Angle {
//	return Angle(NormAngle(a.Radians()))*Radian
//}
