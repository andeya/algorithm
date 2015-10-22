package evaler_test

import (
	"math"
	"math/big"
	"testing"

	"github.com/henrylee2cn/algorithm/evaler"
)

// -----------------------------------------------------------------------------

var testsEval = []struct {
	in  string
	out *big.Rat
	ok  bool
}{
	{"5 + 2", big.NewRat(7, 1), true},            // simple plus
	{"5 - 2", big.NewRat(3, 1), true},            // simple minus
	{"5 * 2", big.NewRat(10, 1), true},           // simple multiply
	{"5 / 2", big.NewRat(5, 2), true},            // simple divide
	{"U + U", nil, false},                        // letters 1
	{"2 + U", nil, false},                        // broken 1
	{"2 +  ", nil, false},                        // broken 2
	{"+ 2 - + * ", nil, false},                   // broken 3
	{"5.5+2*(3+1)", big.NewRat(27, 2), true},     // complex 1
	{"(((1+2.3)))", big.NewRat(33, 10), true},    // complex 2
	{"(1+(2))*(5-2.5)", big.NewRat(15, 2), true}, // complex 3
	{"3*(2<4)", big.NewRat(3, 1), true},          // less than
	{"3*(2>4)", new(big.Rat), true},              // greater than
	{"5 / 0", nil, false},                        // divide by zero
	{"2 ** 3", big.NewRat(8, 1), true},           // exponent 1
	{"9.0**0.5", big.NewRat(3, 1), true},         // exponent 2
	{"1.23", big.NewRat(123, 100), true},
	{"-1+2", big.NewRat(1, 1), true},          // unary minus (the beginning of a expression)
	{"3*-4", big.NewRat(-12, 1), true},        // unary minus (after an operator)
	{"4/(-1+3)", big.NewRat(2, 1), true},      // unary minus (after '(' )
	{"-(-1+2)--2**3", big.NewRat(7, 1), true}, // unary minus (complex)
}

func TestEval(t *testing.T) {
	for i, test := range testsEval {
		ret, err := evaler.Eval(test.in)
		if ret == nil && test.out == nil {
			// ok, do nothing
		} else if ret == nil || test.out == nil {
			t.Errorf("#%d: %s: unexpected nil result: %v vs %v", i, test.in, ret, test.out)
		} else if ret.Cmp(test.out) != 0 {
			t.Errorf("#%d: %s: bad result: got %v expected %v", i, test.in, ret, test.out)
		}
		if (err == nil) != test.ok {
			t.Errorf("#%d: %s: unexpected err result: %t vs %t", i, test.in, (err == nil), test.ok)
		}
	}
}

// -----------------------------------------------------------------------------

var testsBigratToInt = []struct {
	in  *big.Rat
	out int64
	ok  bool
}{
	{big.NewRat(4, 2), int64(2), true},
	{big.NewRat(5, 2), int64(3), true},
	{big.NewRat(-4, 2), int64(-2), true},
	{new(big.Rat).Mul(big.NewRat(math.MaxInt64, 1), big.NewRat(math.MaxInt64, 1)), int64(0), false},
}

func TestBigratToInt(t *testing.T) {
	for i, test := range testsBigratToInt {
		ret, err := evaler.BigratToInt(test.in)
		if test.ok && (ret != test.out) {
			t.Errorf("#%d: got %d expected %d", i, ret, test.out)
		}
		if (err == nil) != test.ok {
			t.Errorf("#%d: %s: unexpected err result: %t vs %t", i, test.in, (err == nil), test.ok)
		}
	}
}

// -----------------------------------------------------------------------------

var testsBigratToBigint = []struct {
	in  *big.Rat
	out *big.Int
}{
	{big.NewRat(4, 2), big.NewInt(2)},
	{big.NewRat(5, 2), big.NewInt(3)},
	{big.NewRat(-4, 2), big.NewInt(-2)},
}

func TestBigratToBigint(t *testing.T) {
	for i, test := range testsBigratToBigint {
		ret := evaler.BigratToBigint(test.in)
		if ret.Cmp(test.out) != 0 {
			t.Errorf("#%d: got %d expected %d", i, ret, test.out)
		}
	}
}

// -----------------------------------------------------------------------------

var testsBigratToFloat = []struct {
	in  *big.Rat
	out float64
}{
	{big.NewRat(4, 2), float64(2.0)},
	{big.NewRat(5, 2), float64(2.5)},
	{big.NewRat(-4, 2), float64(-2.0)},
}

func TestBigratToFloat(t *testing.T) {
	for i, test := range testsBigratToFloat {
		ret := evaler.BigratToFloat(test.in)
		if ret != test.out {
			t.Errorf("#%d: got %d expected %d", i, ret, test.out)
		}
	}
}

// -----------------------------------------------------------------------------

var testsFloatToBigrat = []struct {
	in  float64
	out *big.Rat
}{
	{float64(2.0), big.NewRat(4, 2)},
	{float64(2.5), big.NewRat(5, 2)},
	{float64(-2.0), big.NewRat(-4, 2)},
}

func TestFloatToBigrat(t *testing.T) {
	for i, test := range testsFloatToBigrat {
		ret := evaler.FloatToBigrat(test.in)
		if ret.Cmp(test.out) != 0 {
			t.Errorf("#%d: got %s expected %s", i, ret, test.out)
		}
	}
}
