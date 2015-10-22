package matops

import (
	"github.com/henrylee2cn/algorithm/matrix"
	"math"
	//"fmt"
)

func zeroDims(a *matrix.FloatMatrix) bool {
	return a.Rows() == 0 || a.Cols() == 0
}

/* From LAPACK/dlapy2.f
 *
 * sqrtX2Y2() returns sqrt(x**2+y**2), taking care not to cause unnecessary
 * overflow.
 */
func sqrtX2Y2(x, y float64) float64 {
	xabs := math.Abs(x)
	yabs := math.Abs(y)
	w := xabs
	if yabs > w {
		w = yabs
	}
	z := xabs
	if yabs < z {
		z = yabs
	}
	if z == 0.0 {
		return w
	}
	return w * math.Sqrt(1.0+(z/w)*(z/w))
}

/* From LAPACK/dlarfg.f
 *
 * DLARFG generates a real elementary reflector H of order n, such
 * that
 *
 *       H * ( alpha ) = ( beta ),   H**T * H = I.
 *           (   x   )   (   0  )
 *
 * where alpha and beta are scalars, and x is an (n-1)-element real
 * vector. H is represented in the form
 *
 *       H = I - tau * ( 1 ) * ( 1 v**T ) ,
 *                     ( v )
 *
 * where tau is a real scalar and v is a real (n-1)-element
 * vector.
 *
 * If the elements of x are all zero, then tau = 0 and H is taken to be
 * the unit matrix.
 *
 * Otherwise  1 <= tau <= 2.
 */
func computeHouseholder(a11, x, tau *matrix.FloatMatrix, flags Flags) {

	// norm_x2 = ||x||_2
	norm_x2 := Norm2(x)
	if norm_x2 == 0.0 {
		//a11.SetAt(0, 0, -a11.GetAt(0, 0))
		tau.SetAt(0, 0, 0.0)
		return
	}

	alpha := a11.GetAt(0, 0)
	sign := 1.0
	if math.Signbit(alpha) {
		sign = -1.0
	}
	// beta = -(alpha / |alpha|) * ||alpha x||
	//      = -sign(alpha) * sqrt(alpha**2, norm_x2**2)
	beta := -sign * sqrtX2Y2(alpha, norm_x2)

	// x = x /(a11 - beta)
	InvScale(x, alpha-beta)

	tau.SetAt(0, 0, (beta-alpha)/beta)
	a11.SetAt(0, 0, beta)
}

/* From LAPACK/dlarf.f
 *
 * Applies a real elementary reflector H to a real m by n matrix A,
 * from either the left or the right. H is represented in the form
 *
 *       H = I - tau * ( 1 ) * ( 1 v.T )
 *                     ( v )
 *
 * where tau is a real scalar and v is a real vector.
 *
 * If tau = 0, then H is taken to be the unit matrix.
 *
 * A is /a1\   a1 := a1 - w1
 *      \A2/   A2 := A2 - v*w1
 *             w1 := tau*(a1 + A2.T*v) if side == LEFT
 *                := tau*(a1 + A2*v)   if side == RIGHT
 *
 * Allocates/frees intermediate work space matrix w1.
 */
func applyHouseholder(tau, v, a1, A2 *matrix.FloatMatrix, flags Flags) {

	tval := tau.GetAt(0, 0)
	if tval == 0.0 {
		return
	}
	w1 := a1.Copy()
	if flags&LEFT != 0 {
		// w1 = a1 + A2.T*v
		MVMult(w1, A2, v, 1.0, 1.0, TRANSA)
	} else {
		// w1 = a1 + A2*v
		MVMult(w1, A2, v, 1.0, 1.0, NOTRANS)
	}

	// w1 = tau*w1
	Scale(w1, tval)

	// a1 = a1 - w1
	a1.Minus(w1)

	// A2 = A2 - v*w1
	MVRankUpdate(A2, v, w1, -1.0)
}

/*
 * Applies a real elementary reflector H to a real m by n matrix A,
 * from either the left or the right. H is represented in the form
 *
 *       H = I - tau * ( 1 ) * ( 1 v.T )
 *                     ( v )
 *
 * where tau is a real scalar and v is a real vector.
 *
 * If tau = 0, then H is taken to be the unit matrix.
 *
 * A is /a1\   a1 := a1 - w1
 *      \A2/   A2 := A2 - v*w1
 *             w1 := tau*(a1 + A2.T*v) if side == LEFT
 *                := tau*(a1 + A2*v)   if side == RIGHT
 *
 * Intermediate work space w1 required as parameter, no allocation.
 */
func applyHHTo2x1(tau, v, a1, A2, w1 *matrix.FloatMatrix, flags Flags) {

	tval := tau.GetAt(0, 0)
	if tval == 0.0 {
		return
	}

	// maybe with Scale(0.0), Axpy(w1, a1, 1.0)
	a1.CopyTo(w1)
	if flags&LEFT != 0 {
		// w1 = a1 + A2.T*v
		MVMult(w1, A2, v, 1.0, 1.0, TRANSA)
	} else {
		// w1 = a1 + A2*v
		MVMult(w1, A2, v, 1.0, 1.0, NOTRANS)
	}

	// w1 = tau*w1
	Scale(w1, tval)

	// a1 = a1 - w1
	a1.Minus(w1)

	// A2 = A2 - v*w1
	if flags&LEFT != 0 {
		MVRankUpdate(A2, v, w1, -1.0)
	} else {
		MVRankUpdate(A2, w1, v, -1.0)
	}
}

func applyHHTo1x1(tau, v, A2, w1 *matrix.FloatMatrix, flags Flags) {

	tval := tau.GetAt(0, 0)
	if tval == 0.0 {
		return
	}
	if flags&LEFT != 0 {
		// w1 = A2.T*v
		MVMult(w1, A2, v, 1.0, 0.0, TRANSA)
	} else {
		// w1 = A2*v
		MVMult(w1, A2, v, 1.0, 0.0, NOTRANS)
	}

	// A2 = A2 - tau*v*w1
	MVRankUpdate(A2, v, w1, -tval)
}

// Local Variables:
// tab-width: 4
// indent-tabs-mode: nil
// End:
