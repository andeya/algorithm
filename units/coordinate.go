package units

import (
	"math"
)

type Coordinate struct {
	X, Y Distance
}

func CoordinateAdd(c1, c2 Coordinate) Coordinate {
	dx := c1.X + c2.X
	dy := c1.Y + c2.Y
	return Coordinate{dx, dy}
}
func CoordinateSubtract(c1, c2 Coordinate) Coordinate {
	dx := c1.X - c2.X
	dy := c1.Y - c2.Y
	return Coordinate{dx, dy}
}

func (c1 Coordinate) Distance(c2 Coordinate) Distance {
	dx := c2.X - c1.X
	dy := c2.Y - c1.Y
	return Distance(math.Hypot(float64(dx), float64(dy)))
}

func (c1 Coordinate) Angle(c2 Coordinate) Angle {
	dx := c2.X - c1.X
	dy := c2.Y - c1.Y
	return Atan2(dy, dx).Normalize()
}
