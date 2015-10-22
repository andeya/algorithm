package fixed

import (
	"math"
	"strconv"
)

const (
	fracBits  = 32
	scale     = 1 << fracBits
	halfScale = scale >> 1
	lowMask   = scale - 1
)

var (
	Zero = Fixed{0}
	One  = Fixed{scale}
	Min  = Fixed{math.MinInt64}
	Max  = Fixed{math.MaxInt64}
)

// Fixed-point number in Q32.32 format.
type Fixed struct {
	value int64
}

func New(x float64) Fixed {
	return Fixed{int64(x * scale)}
}

func (x Fixed) Abs() Fixed {
	if x.value < 0 {
		return x.Neg()
	}
	return x
}

func (x Fixed) Neg() Fixed {
	return Fixed{-x.value}
}

func (x Fixed) Add(y Fixed) Fixed {
	return Fixed{x.value + y.value}
}

func (x Fixed) Addf(y float64) Fixed {
	return x.Add(New(y))
}

func (x Fixed) Sub(y Fixed) Fixed {
	return Fixed{x.value - y.value}
}

func (x Fixed) Subf(y float64) Fixed {
	return x.Sub(New(y))
}

func (x Fixed) Mul(y Fixed) Fixed {
	return Fixed{(x.value / halfScale) * (y.value / halfScale)}
}

func (x Fixed) Mulf(y float64) Fixed {
	return x.Mul(New(y))
}

func (x Fixed) Div(y Fixed) Fixed {
	return Fixed{(x.value * halfScale) / (y.value * halfScale)}
}

func (x Fixed) Divf(y float64) Fixed {
	return x.Div(New(y))
}

func (x Fixed) Mod(y Fixed) Fixed {
	return Fixed{x.value % y.value}
}

func (x Fixed) Modf(y float64) Fixed {
	return x.Mod(New(y))
}

func (x Fixed) Ceil() Fixed {
	if x.value&lowMask == 0 {
		return x
	}
	return Fixed{x.value&^lowMask + scale}
}

func (x Fixed) Floor() Fixed {
	return Fixed{x.value &^ lowMask}
}

func (x Fixed) Gt(y Fixed) bool {
	return x.value > y.value
}

func (x Fixed) Lt(y Fixed) bool {
	return x.value < y.value
}

func (x Fixed) Geq(y Fixed) bool {
	return x.value >= y.value
}

func (x Fixed) Leq(y Fixed) bool {
	return x.value <= y.value
}

func (x Fixed) Cmp(y Fixed) int {
	if x.value > y.value {
		return 1
	}
	if x.value < y.value {
		return -1
	}
	return 0
}

func (x Fixed) Max(y Fixed) Fixed {
	if x.value >= y.value {
		return x
	}
	return y
}

func (x Fixed) Min(y Fixed) Fixed {
	if x.value <= y.value {
		return x
	}
	return y
}

func (x Fixed) Int64() int64 {
	if x.value >= 0 || x.value&lowMask == 0 {
		return int64(x.value >> fracBits)
	}
	return int64(x.value>>fracBits) + 1
}

func (x Fixed) Int32() int32 {
	return int32(x.Int64())
}

func (x Fixed) Frac() float64 {
	return float64(x.value & lowMask)
}

func (x Fixed) Float64() float64 {
	return float64(x.value) / scale
}

func (x Fixed) String() string {
	// XXX This is an ugly hack.
	int := x.Int64()
	frac := x.Frac()
	if x.value >= 0 {
		return strconv.FormatInt(int, 10) + strconv.FormatFloat(frac, 'f', -1, 64)[1:]
	}
	if int == 0 {
		return "-" + strconv.FormatInt(int, 10) + strconv.FormatFloat(frac, 'f', -1, 64)[2:]
	}
	return strconv.FormatInt(int, 10) + strconv.FormatFloat(frac, 'f', -1, 64)[2:]
}
