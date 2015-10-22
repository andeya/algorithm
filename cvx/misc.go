// Copyright (c) Harri Rautila, 2012

// This file is part of github.com/hrautila/cvx package.
// It is free software, distributed under the terms of GNU Lesser General Public
// License Version 3, or any later version. See the COPYING tile included in this archive.

package cvx

import (
	"github.com/henrylee2cn/algorithm/cvx/sets"
	la_ "github.com/henrylee2cn/algorithm/linalg"
	"github.com/henrylee2cn/algorithm/linalg/blas"
	"github.com/henrylee2cn/algorithm/linalg/lapack"
	"github.com/henrylee2cn/algorithm/matrix"
	//"github.com/henrylee2cn/algorithm/cvx/checkpnt"
	"math"
	//"fmt"
)

func maxdim(vec []int) int {
	res := 0
	for _, v := range vec {
		if v > res {
			res = v
		}
	}
	return res
}

func sumdim(vec []int) int {
	res := 0
	for _, v := range vec {
		res += v
	}
	return res
}

func maxvec(vec []float64) float64 {
	res := math.Inf(-1)
	for _, v := range vec {
		if v > res {
			res = v
		}
	}
	return res
}

func minvec(vec []float64) float64 {
	res := math.Inf(+1)
	for _, v := range vec {
		if v < res {
			res = v
		}
	}
	return res
}

/*
   Applies Nesterov-Todd scaling or its inverse.

   Computes

        x := W*x        (trans is false 'N', inverse = false 'N')
        x := W^T*x      (trans is true  'T', inverse = false 'N')
        x := W^{-1}*x   (trans is false 'N', inverse = true  'T')
        x := W^{-T}*x   (trans is true  'T', inverse = true  'T').

   x is a dense float matrix.

   W is a MatrixSet with entries:

   - W['dnl']: positive vector
   - W['dnli']: componentwise inverse of W['dnl']
   - W['d']: positive vector
   - W['di']: componentwise inverse of W['d']
   - W['v']: lists of 2nd order cone vectors with unit hyperbolic norms
   - W['beta']: list of positive numbers
   - W['r']: list of square matrices
   - W['rti']: list of square matrices.  rti[k] is the inverse transpose
     of r[k].

   The 'dnl' and 'dnli' entries are optional, and only present when the
   function is called from the nonlinear solver.
*/
func scale(x *matrix.FloatMatrix, W *sets.FloatMatrixSet, trans, inverse bool) (err error) {
	/*DEBUGGED*/
	var wl []*matrix.FloatMatrix
	var w *matrix.FloatMatrix
	ind := 0
	err = nil

	// var minor int = 0
	//if ! checkpnt.MinorEmpty() {
	//	minor = checkpnt.MinorTop()
	//}

	//fmt.Printf("\n%d.%04d scaling x=\n%v\n", checkpnt.Major(), minor, x.ToString("%.17f"))

	// Scaling for nonlinear component xk is xk := dnl .* xk; inverse
	// scaling is xk ./ dnl = dnli .* xk, where dnl = W['dnl'],
	// dnli = W['dnli'].

	if wl = W.At("dnl"); wl != nil {
		if inverse {
			w = W.At("dnli")[0]
		} else {
			w = W.At("dnl")[0]
		}
		for k := 0; k < x.Cols(); k++ {
			err = blas.TbmvFloat(w, x, &la_.IOpt{"n", w.Rows()}, &la_.IOpt{"k", 0},
				&la_.IOpt{"lda", 1}, &la_.IOpt{"offsetx", k * x.Rows()})
			if err != nil {
				//fmt.Printf("1. TbmvFloat: %v\n", err)
				return
			}
		}
		ind += w.Rows()
	}

	//if ! checkpnt.MinorEmpty() {
	//    checkpnt.Check("000scale", minor)
	//}

	// Scaling for linear 'l' component xk is xk := d .* xk; inverse
	// scaling is xk ./ d = di .* xk, where d = W['d'], di = W['di'].

	if inverse {
		w = W.At("di")[0]
	} else {
		w = W.At("d")[0]
	}

	for k := 0; k < x.Cols(); k++ {
		err = blas.TbmvFloat(w, x, &la_.IOpt{"n", w.Rows()}, &la_.IOpt{"k", 0},
			&la_.IOpt{"lda", 1}, &la_.IOpt{"offsetx", k*x.Rows() + ind})
		if err != nil {
			//fmt.Printf("2. TbmvFloat: %v\n", err)
			return
		}
	}
	ind += w.Rows()

	//if ! checkpnt.MinorEmpty() {
	//	checkpnt.Check("010scale", minor)
	//}

	// Scaling for 'q' component is
	//
	//    xk := beta * (2*v*v' - J) * xk
	//        = beta * (2*v*(xk'*v)' - J*xk)
	//
	// where beta = W['beta'][k], v = W['v'][k], J = [1, 0; 0, -I].
	//
	//Inverse scaling is
	//
	//    xk := 1/beta * (2*J*v*v'*J - J) * xk
	//        = 1/beta * (-J) * (2*v*((-J*xk)'*v)' + xk).
	//wf := matrix.FloatZeros(x.Cols(), 1)
	w = matrix.FloatZeros(x.Cols(), 1)
	for k, v := range W.At("v") {
		m := v.Rows()
		if inverse {
			blas.ScalFloat(x, -1.0, &la_.IOpt{"offset", ind}, &la_.IOpt{"inc", x.Rows()})
		}
		err = blas.GemvFloat(x, v, w, 1.0, 0.0, la_.OptTrans, &la_.IOpt{"m", m},
			&la_.IOpt{"n", x.Cols()}, &la_.IOpt{"offsetA", ind},
			&la_.IOpt{"lda", x.Rows()})
		if err != nil {
			//fmt.Printf("3. GemvFloat: %v\n", err)
			return
		}

		err = blas.ScalFloat(x, -1.0, &la_.IOpt{"offset", ind}, &la_.IOpt{"inc", x.Rows()})
		if err != nil {
			return
		}

		err = blas.GerFloat(v, w, x, 2.0, &la_.IOpt{"m", m},
			&la_.IOpt{"n", x.Cols()}, &la_.IOpt{"lda", x.Rows()},
			&la_.IOpt{"offsetA", ind})
		if err != nil {
			//fmt.Printf("4. GerFloat: %v\n", err)
			return
		}

		var a float64
		if inverse {
			blas.ScalFloat(x, -1.0,
				&la_.IOpt{"offset", ind}, &la_.IOpt{"inc", x.Rows()})
			// a[i,j] := 1.0/W[i,j]
			a = 1.0 / W.At("beta")[0].GetIndex(k)
		} else {
			a = W.At("beta")[0].GetIndex(k)
		}
		for i := 0; i < x.Cols(); i++ {
			blas.ScalFloat(x, a, &la_.IOpt{"n", m}, &la_.IOpt{"offset", ind + i*x.Rows()})
		}
		ind += m
	}

	//if ! checkpnt.MinorEmpty() {
	//	checkpnt.Check("020scale", minor)
	//}

	// Scaling for 's' component xk is
	//
	//     xk := vec( r' * mat(xk) * r )  if trans = 'N'
	//     xk := vec( r * mat(xk) * r' )  if trans = 'T'.
	//
	// r is kth element of W['r'].
	//
	// Inverse scaling is
	//
	//     xk := vec( rti * mat(xk) * rti' )  if trans = 'N'
	//     xk := vec( rti' * mat(xk) * rti )  if trans = 'T'.
	//
	// rti is kth element of W['rti'].
	maxn := 0
	for _, r := range W.At("r") {
		if r.Rows() > maxn {
			maxn = r.Rows()
		}
	}
	a := matrix.FloatZeros(maxn, maxn)
	for k, v := range W.At("r") {
		t := trans
		var r *matrix.FloatMatrix
		if !inverse {
			r = v
			t = !trans
		} else {
			r = W.At("rti")[k]
		}

		n := r.Rows()
		for i := 0; i < x.Cols(); i++ {
			// scale diagonal of xk by 0.5
			blas.ScalFloat(x, 0.5, &la_.IOpt{"offset", ind + i*x.Rows()},
				&la_.IOpt{"inc", n + 1}, &la_.IOpt{"n", n})

			// a = r*tril(x) (t is 'N') or a = tril(x)*r  (t is 'T')
			blas.Copy(r, a)
			if !t {
				err = blas.TrmmFloat(x, a, 1.0, la_.OptRight, &la_.IOpt{"m", n},
					&la_.IOpt{"n", n}, &la_.IOpt{"lda", n}, &la_.IOpt{"ldb", n},
					&la_.IOpt{"offsetA", ind + i*x.Rows()})
				if err != nil {
					//fmt.Printf("5. TrmmFloat: %v\n", err)
					return
				}

				// x := (r*a' + a*r')  if t is 'N'
				err = blas.Syr2kFloat(r, a, x, 1.0, 0.0, la_.OptNoTrans, &la_.IOpt{"n", n},
					&la_.IOpt{"k", n}, &la_.IOpt{"ldb", n}, &la_.IOpt{"ldc", n},
					&la_.IOpt{"offsetC", ind + i*x.Rows()})
				if err != nil {
					//fmt.Printf("6. Syr2kFloat: %v\n", err)
					return
				}

			} else {
				err = blas.TrmmFloat(x, a, 1.0, la_.OptLeft, &la_.IOpt{"m", n},
					&la_.IOpt{"n", n}, &la_.IOpt{"lda", n}, &la_.IOpt{"ldb", n},
					&la_.IOpt{"offsetA", ind + i*x.Rows()})
				if err != nil {
					//fmt.Printf("7. TrmmFloat: %v\n", err)
					return
				}

				// x := (r'*a + a'*r)  if t is 'T'
				err = blas.Syr2kFloat(r, a, x, 1.0, 0.0, la_.OptTrans, &la_.IOpt{"n", n},
					&la_.IOpt{"k", n}, &la_.IOpt{"ldb", n}, &la_.IOpt{"ldc", n},
					&la_.IOpt{"offsetC", ind + i*x.Rows()})
				if err != nil {
					//fmt.Printf("8. Syr2kFloat: %v\n", err)
					return
				}
			}
		}
		ind += n * n
	}
	//if ! checkpnt.MinorEmpty() {
	//	checkpnt.Check("030scale", minor)
	//}
	return
}

/*
   Evaluates

       x := H(lambda^{1/2}) * x   (inverse is 'N')
       x := H(lambda^{-1/2}) * x  (inverse is 'I').

   H is the Hessian of the logarithmic barrier.

*/
func scale2(lmbda, x *matrix.FloatMatrix, dims *sets.DimensionSet, mnl int, inverse bool) (err error) {
	err = nil

	//var minor int = 0
	//if ! checkpnt.MinorEmpty() {
	//	minor = checkpnt.MinorTop()
	//}

	//fmt.Printf("\n%d.%04d scale2 x=\n%v\nlmbda=\n%v\n", checkpnt.Major(), minor,
	//	x.ToString("%.17f"), lmbda.ToString("%.17f"))

	//if ! checkpnt.MinorEmpty() {
	//	checkpnt.Check("000scale2", minor)
	//}

	// For the nonlinear and 'l' blocks,
	//
	//     xk := xk ./ l   (inverse is 'N')
	//     xk := xk .* l   (inverse is 'I')
	//
	// where l is lmbda[:mnl+dims['l']].
	ind := mnl + dims.Sum("l")
	if !inverse {
		blas.TbsvFloat(lmbda, x, &la_.IOpt{"n", ind}, &la_.IOpt{"k", 0}, &la_.IOpt{"lda", 1})
	} else {
		blas.TbmvFloat(lmbda, x, &la_.IOpt{"n", ind}, &la_.IOpt{"k", 0}, &la_.IOpt{"lda", 1})
	}

	//if ! checkpnt.MinorEmpty() {
	//	checkpnt.Check("010scale2", minor)
	//}

	// For 'q' blocks, if inverse is 'N',
	//
	//     xk := 1/a * [ l'*J*xk;
	//         xk[1:] - (xk[0] + l'*J*xk) / (l[0] + 1) * l[1:] ].
	//
	// If inverse is 'I',
	//
	//     xk := a * [ l'*xk;
	//         xk[1:] + (xk[0] + l'*xk) / (l[0] + 1) * l[1:] ].
	//
	// a = sqrt(lambda_k' * J * lambda_k), l = lambda_k / a.
	for _, m := range dims.At("q") {
		var lx, a, c, x0 float64
		a = jnrm2(lmbda, m, ind) //&la_.IOpt{"n", m}, &la_.IOpt{"offset", ind})
		if !inverse {
			lx = jdot(lmbda, x, m, ind, ind) //&la_.IOpt{"n", m}, &la_.IOpt{"offsetx", ind},
			//&la_.IOpt{"offsety", ind})
			lx /= a
		} else {
			lx = blas.DotFloat(lmbda, x, &la_.IOpt{"n", m}, &la_.IOpt{"offsetx", ind},
				&la_.IOpt{"offsety", ind})
			lx /= a
		}
		x0 = x.GetIndex(ind)
		x.SetIndex(ind, lx)
		c = (lx + x0) / (lmbda.GetIndex(ind)/a + 1.0) / a
		if !inverse {
			c *= -1.0
		}
		blas.AxpyFloat(lmbda, x, c, &la_.IOpt{"n", m - 1}, &la_.IOpt{"offsetx", ind + 1},
			&la_.IOpt{"offsety", ind + 1})
		if !inverse {
			a = 1.0 / a
		}
		blas.ScalFloat(x, a, &la_.IOpt{"offset", ind}, &la_.IOpt{"n", m})
		ind += m
	}

	//if ! checkpnt.MinorEmpty() {
	//	checkpnt.Check("020scale2", minor)
	//}

	// For the 's' blocks, if inverse is 'N',
	//
	//     xk := vec( diag(l)^{-1/2} * mat(xk) * diag(k)^{-1/2}).
	//
	// If inverse is true,
	//
	//     xk := vec( diag(l)^{1/2} * mat(xk) * diag(k)^{1/2}).
	//
	// where l is kth block of lambda.
	//
	// We scale upper and lower triangular part of mat(xk) because the
	// inverse operation will be applied to nonsymmetric matrices.
	ind2 := ind
	sdims := dims.At("s")
	for k := 0; k < len(sdims); k++ {
		m := sdims[k]
		scaleF := func(v, x float64) float64 {
			return math.Sqrt(v) * math.Sqrt(x)
		}
		for j := 0; j < m; j++ {
			c := matrix.FloatVector(lmbda.FloatArray()[ind2 : ind2+m])
			c.ApplyConst(lmbda.GetIndex(ind2+j), scaleF)
			if !inverse {
				blas.Tbsv(c, x, &la_.IOpt{"n", m}, &la_.IOpt{"k", 0}, &la_.IOpt{"lda", 1},
					&la_.IOpt{"offsetx", ind + j*m})
			} else {
				blas.Tbmv(c, x, &la_.IOpt{"n", m}, &la_.IOpt{"k", 0}, &la_.IOpt{"lda", 1},
					&la_.IOpt{"offsetx", ind + j*m})
			}
		}
		ind += m * m
		ind2 += m
	}

	//if ! checkpnt.MinorEmpty() {
	//	checkpnt.Check("030scale2", minor)
	//}
	return
}

/*
   Updates the Nesterov-Todd scaling matrix W and the scaled variable
   lmbda so that on exit

         W * zt = W^{-T} * st = lmbda.

   On entry, the nonlinear, 'l' and 'q' components of the arguments s
   and z contain W^{-T}*st and W*zt, i.e, the new iterates in the current
   scaling.

   The 's' components contain the factors Ls, Lz in a factorization of
   the new iterates in the current scaling, W^{-T}*st = Ls*Ls',
   W*zt = Lz*Lz'.

*/

func updateScaling(W *sets.FloatMatrixSet, lmbda, s, z *matrix.FloatMatrix) (err error) {
	err = nil
	var stmp, ztmp *matrix.FloatMatrix
	/*
	   Nonlinear and 'l' blocks

	      d :=  d .* sqrt( s ./ z )
	      lmbda := lmbda .* sqrt(s) .* sqrt(z)
	*/
	mnl := 0
	dnlset := W.At("dnl")
	dnliset := W.At("dnli")
	dset := W.At("d")
	diset := W.At("di")
	beta := W.At("beta")[0]
	if dnlset != nil && dnlset[0].NumElements() > 0 {
		mnl = dnlset[0].NumElements()
	}
	ml := dset[0].NumElements()
	m := mnl + ml
	//fmt.Printf("ml=%d, mnl=%d, m=%d'n", ml, mnl, m)

	stmp = matrix.FloatVector(s.FloatArray()[:m])
	stmp.Apply(math.Sqrt)
	s.SetIndexesFromArray(stmp.FloatArray(), matrix.MakeIndexSet(0, m, 1)...)

	ztmp = matrix.FloatVector(z.FloatArray()[:m])
	ztmp.Apply(math.Sqrt)
	z.SetIndexesFromArray(ztmp.FloatArray(), matrix.MakeIndexSet(0, m, 1)...)

	// d := d .* s .* z
	if len(dnlset) > 0 {
		blas.TbmvFloat(s, dnlset[0], &la_.IOpt{"n", mnl}, &la_.IOpt{"k", 0}, &la_.IOpt{"lda", 1})
		blas.TbsvFloat(z, dnlset[0], &la_.IOpt{"n", mnl}, &la_.IOpt{"k", 0}, &la_.IOpt{"lda", 1})
		//dnliset[0].Apply(dnlset[0], func(a float64)float64 { return 1.0/a})
		//--dnliset[0] = matrix.Inv(dnlset[0])
		matrix.Set(dnliset[0], dnlset[0])
		dnliset[0].Inv()
	}
	blas.TbmvFloat(s, dset[0], &la_.IOpt{"n", ml},
		&la_.IOpt{"k", 0}, &la_.IOpt{"lda", 1}, &la_.IOpt{"offseta", mnl})
	blas.TbsvFloat(z, dset[0], &la_.IOpt{"n", ml},
		&la_.IOpt{"k", 0}, &la_.IOpt{"lda", 1}, &la_.IOpt{"offseta", mnl})
	//diset[0].Apply(dset[0], func(a float64)float64 { return 1.0/a})
	//--diset[0] = matrix.Inv(dset[0])
	matrix.Set(diset[0], dset[0])
	diset[0].Inv()

	// lmbda := s .* z
	blas.CopyFloat(s, lmbda, &la_.IOpt{"n", m})
	blas.TbmvFloat(z, lmbda, &la_.IOpt{"n", m}, &la_.IOpt{"k", 0}, &la_.IOpt{"lda", 1})

	// 'q' blocks.
	// Let st and zt be the new variables in the old scaling:
	//
	//     st = s_k,   zt = z_k
	//
	// and a = sqrt(st' * J * st),  b = sqrt(zt' * J * zt).
	//
	// 1. Compute the hyperbolic Householder transformation 2*q*q' - J
	//    that maps st/a to zt/b.
	//
	//        c = sqrt( (1 + st'*zt/(a*b)) / 2 )
	//        q = (st/a + J*zt/b) / (2*c).
	//
	//    The new scaling point is
	//
	//        wk := betak * sqrt(a/b) * (2*v[k]*v[k]' - J) * q
	//
	//    with betak = W['beta'][k].
	//
	// 3. The scaled variable:
	//
	//        lambda_k0 = sqrt(a*b) * c
	//        lambda_k1 = sqrt(a*b) * ( (2vk*vk' - J) * (-d*q + u/2) )_1
	//
	//    where
	//
	//        u = st/a - J*zt/b
	//        d = ( vk0 * (vk'*u) + u0/2 ) / (2*vk0 *(vk'*q) - q0 + 1).
	//
	// 4. Update scaling
	//
	//        v[k] := wk^1/2
	//              = 1 / sqrt(2*(wk0 + 1)) * (wk + e).
	//        beta[k] *=  sqrt(a/b)

	ind := m
	for k, v := range W.At("v") {
		m = v.NumElements()

		// ln = sqrt( lambda_k' * J * lambda_k ) !! NOT USED!!
		jnrm2(lmbda, m, ind) // ?? NOT USED ??

		// a = sqrt( sk' * J * sk ) = sqrt( st' * J * st )
		// s := s / a = st / a
		aa := jnrm2(s, m, ind)
		blas.ScalFloat(s, 1.0/aa, &la_.IOpt{"n", m}, &la_.IOpt{"offset", ind})

		// b = sqrt( zk' * J * zk ) = sqrt( zt' * J * zt )
		// z := z / a = zt / b
		bb := jnrm2(z, m, ind)
		blas.ScalFloat(z, 1.0/bb, &la_.IOpt{"n", m}, &la_.IOpt{"offset", ind})

		// c = sqrt( ( 1 + (st'*zt) / (a*b) ) / 2 )
		cc := blas.DotFloat(s, z, &la_.IOpt{"offsetx", ind}, &la_.IOpt{"offsety", ind},
			&la_.IOpt{"n", m})
		cc = math.Sqrt((1.0 + cc) / 2.0)

		// vs = v' * st / a
		vs := blas.DotFloat(v, s, &la_.IOpt{"offsety", ind}, &la_.IOpt{"n", m})

		// vz = v' * J *zt / b
		vz := jdot(v, z, m, 0, ind)

		// vq = v' * q where q = (st/a + J * zt/b) / (2 * c)
		vq := (vs + vz) / 2.0 / cc

		// vq = v' * q where q = (st/a + J * zt/b) / (2 * c)
		vu := vs - vz
		// lambda_k0 = c
		lmbda.SetIndex(ind, cc)

		// wk0 = 2 * vk0 * (vk' * q) - q0
		wk0 := 2.0*v.GetIndex(0)*vq - (s.GetIndex(ind)+z.GetIndex(ind))/2.0/cc

		// d = (v[0] * (vk' * u) - u0/2) / (wk0 + 1)
		dd := (v.GetIndex(0)*vu - s.GetIndex(ind)/2.0 + z.GetIndex(ind)/2.0) / (wk0 + 1.0)

		// lambda_k1 = 2 * v_k1 * vk' * (-d*q + u/2) - d*q1 + u1/2
		blas.CopyFloat(v, lmbda, &la_.IOpt{"offsetx", 1}, &la_.IOpt{"offsety", ind + 1},
			&la_.IOpt{"n", m - 1})
		blas.ScalFloat(lmbda, (2.0 * (-dd*vq + 0.5*vu)),
			&la_.IOpt{"offsetx", ind + 1}, &la_.IOpt{"offsety", ind + 1}, &la_.IOpt{"n", m - 1})
		blas.AxpyFloat(s, lmbda, 0.5*(1.0-dd/cc),
			&la_.IOpt{"offsetx", ind + 1}, &la_.IOpt{"offsety", ind + 1}, &la_.IOpt{"n", m - 1})
		blas.AxpyFloat(z, lmbda, 0.5*(1.0+dd/cc),
			&la_.IOpt{"offsetx", ind + 1}, &la_.IOpt{"offsety", ind + 1}, &la_.IOpt{"n", m - 1})

		// Scale so that sqrt(lambda_k' * J * lambda_k) = sqrt(aa*bb).
		blas.ScalFloat(lmbda, math.Sqrt(aa*bb), &la_.IOpt{"offset", ind}, &la_.IOpt{"n", m})

		// v := (2*v*v' - J) * q
		//    = 2 * (v'*q) * v' - (J* st/a + zt/b) / (2*c)
		blas.ScalFloat(v, 2.0*vq)
		v.SetIndex(0, v.GetIndex(0)-(s.GetIndex(ind)/2.0/cc))
		blas.AxpyFloat(s, v, 0.5/cc, &la_.IOpt{"offsetx", ind + 1}, &la_.IOpt{"offsety", 1},
			&la_.IOpt{"n", m - 1})
		blas.AxpyFloat(z, v, -0.5/cc, &la_.IOpt{"offsetx", ind}, &la_.IOpt{"n", m})

		// v := v^{1/2} = 1/sqrt(2 * (v0 + 1)) * (v + e)
		v0 := v.GetIndex(0) + 1.0
		v.SetIndex(0, v0)
		blas.ScalFloat(v, 1.0/math.Sqrt(2.0*v0))

		// beta[k] *= ( aa / bb )**1/2
		bk := beta.GetIndex(k)
		beta.SetIndex(k, bk*math.Sqrt(aa/bb))

		ind += m
	}
	//fmt.Printf("-- end of q:\nz=\n%v\nlmbda=\n%v\n", z.ConvertToString(), lmbda.ConvertToString())
	//fmt.Printf("beta=\n%v\n", beta.ConvertToString())

	// 's' blocks
	//
	// Let st, zt be the updated variables in the old scaling:
	//
	//     st = Ls * Ls', zt = Lz * Lz'.
	//
	// where Ls and Lz are the 's' components of s, z.
	//
	// 1.  SVD Lz'*Ls = Uk * lambda_k^+ * Vk'.
	//
	// 2.  New scaling is
	//
	//         r[k] := r[k] * Ls * Vk * diag(lambda_k^+)^{-1/2}
	//         rti[k] := r[k] * Lz * Uk * diag(lambda_k^+)^{-1/2}.
	//

	maxr := 0
	for _, m := range W.At("r") {
		if m.Rows() > maxr {
			maxr = m.Rows()
		}
	}
	work := matrix.FloatZeros(maxr*maxr, 1)
	vlensum := 0
	for _, m := range W.At("v") {
		vlensum += m.NumElements()
	}
	ind = mnl + ml + vlensum
	ind2 := ind
	ind3 := 0
	rset := W.At("r")
	rtiset := W.At("rti")

	for k, _ := range rset {
		r := rset[k]
		rti := rtiset[k]
		m = r.Rows()
		//fmt.Printf("m=%d, r=\n%v\nrti=\n%v\n", m, r.ConvertToString(), rti.ConvertToString())

		// r := r*sk = r*Ls
		blas.GemmFloat(r, s, work, 1.0, 0.0, &la_.IOpt{"m", m}, &la_.IOpt{"n", m},
			&la_.IOpt{"k", m}, &la_.IOpt{"ldb", m}, &la_.IOpt{"ldc", m},
			&la_.IOpt{"offsetb", ind2})
		//fmt.Printf("1 work=\n%v\n", work.ConvertToString())
		blas.CopyFloat(work, r, &la_.IOpt{"n", m * m})

		// rti := rti*zk = rti*Lz
		blas.GemmFloat(rti, z, work, 1.0, 0.0, &la_.IOpt{"m", m}, &la_.IOpt{"n", m},
			&la_.IOpt{"k", m}, &la_.IOpt{"ldb", m}, &la_.IOpt{"ldc", m},
			&la_.IOpt{"offsetb", ind2})
		//fmt.Printf("2 work=\n%v\n", work.ConvertToString())
		blas.CopyFloat(work, rti, &la_.IOpt{"n", m * m})

		// SVD Lz'*Ls = U * lmbds^+ * V'; store U in sk and V' in zk. '
		blas.GemmFloat(z, s, work, 1.0, 0.0, la_.OptTransA, &la_.IOpt{"m", m},
			&la_.IOpt{"n", m}, &la_.IOpt{"k", m}, &la_.IOpt{"lda", m}, &la_.IOpt{"ldb", m},
			&la_.IOpt{"ldc", m}, &la_.IOpt{"offseta", ind2}, &la_.IOpt{"offsetb", ind2})
		//fmt.Printf("3 work=\n%v\n", work.ConvertToString())

		// U = s, Vt = z
		lapack.GesvdFloat(work, lmbda, s, z, la_.OptJobuAll, la_.OptJobvtAll,
			&la_.IOpt{"m", m}, &la_.IOpt{"n", m}, &la_.IOpt{"lda", m}, &la_.IOpt{"ldu", m},
			&la_.IOpt{"ldvt", m}, &la_.IOpt{"offsets", ind}, &la_.IOpt{"offsetu", ind2},
			&la_.IOpt{"offsetvt", ind2})

		// r := r*V
		blas.GemmFloat(r, z, work, 1.0, 0.0, la_.OptTransB, &la_.IOpt{"m", m},
			&la_.IOpt{"n", m}, &la_.IOpt{"k", m}, &la_.IOpt{"ldb", m}, &la_.IOpt{"ldc", m},
			&la_.IOpt{"offsetb", ind2})
		//fmt.Printf("4 work=\n%v\n", work.ConvertToString())
		blas.CopyFloat(work, r, &la_.IOpt{"n", m * m})

		// rti := rti*U
		blas.GemmFloat(rti, s, work, 1.0, 0.0, &la_.IOpt{"m", m}, &la_.IOpt{"n", m},
			&la_.IOpt{"k", m}, &la_.IOpt{"ldb", m}, &la_.IOpt{"ldc", m},
			&la_.IOpt{"offsetb", ind2})
		//fmt.Printf("5 work=\n%v\n", work.ConvertToString())
		blas.CopyFloat(work, rti, &la_.IOpt{"n", m * m})

		for i := 0; i < m; i++ {
			a := 1.0 / math.Sqrt(lmbda.GetIndex(ind+i))
			blas.ScalFloat(r, a, &la_.IOpt{"n", m}, &la_.IOpt{"offset", m * i})
			blas.ScalFloat(rti, a, &la_.IOpt{"n", m}, &la_.IOpt{"offset", m * i})
		}
		ind += m
		ind2 += m * m
		ind3 += m // !!NOT USED: ind3!!
	}

	//fmt.Printf("-- end of s:\nz=\n%v\nlmbda=\n%v\n", z.ConvertToString(), lmbda.ConvertToString())

	return

}

/*
   Returns the Nesterov-Todd scaling W at points s and z, and stores the
   scaled variable in lmbda.

       W * z = W^{-T} * s = lmbda.

   W is a MatrixSet with entries:

   - W['dnl']: positive vector
   - W['dnli']: componentwise inverse of W['dnl']
   - W['d']: positive vector
   - W['di']: componentwise inverse of W['d']
   - W['v']: lists of 2nd order cone vectors with unit hyperbolic norms
   - W['beta']: list of positive numbers
   - W['r']: list of square matrices
   - W['rti']: list of square matrices.  rti[k] is the inverse transpose
     of r[k].

*/
func computeScaling(s, z, lmbda *matrix.FloatMatrix, dims *sets.DimensionSet, mnl int) (W *sets.FloatMatrixSet, err error) {
	/*DEBUGGED*/
	err = nil
	W = sets.NewFloatSet("dnl", "dnli", "d", "di", "v", "beta", "r", "rti")

	// For the nonlinear block:
	//
	//     W['dnl'] = sqrt( s[:mnl] ./ z[:mnl] )
	//     W['dnli'] = sqrt( z[:mnl] ./ s[:mnl] )
	//     lambda[:mnl] = sqrt( s[:mnl] .* z[:mnl] )

	var stmp, ztmp, lmd *matrix.FloatMatrix
	if mnl > 0 {
		stmp = matrix.FloatVector(s.FloatArray()[:mnl])
		ztmp = matrix.FloatVector(z.FloatArray()[:mnl])
		//dnl := stmp.Div(ztmp)
		//dnl.Apply(dnl, math.Sqrt)
		dnl := matrix.Sqrt(matrix.Div(stmp, ztmp))
		//dnli := dnl.Copy()
		//dnli.Apply(dnli, func(a float64)float64 { return 1.0/a })
		dnli := matrix.Inv(dnl)
		W.Set("dnl", dnl)
		W.Set("dnli", dnli)
		//lmd = stmp.Mul(ztmp)
		//lmd.Apply(lmd, math.Sqrt)
		lmd = matrix.Sqrt(matrix.Mul(stmp, ztmp))
		lmbda.SetIndexesFromArray(lmd.FloatArray(), matrix.MakeIndexSet(0, mnl, 1)...)
	} else {
		// set for empty matrices
		//W.Set("dnl", matrix.FloatZeros(0, 1))
		//W.Set("dnli", matrix.FloatZeros(0, 1))
		mnl = 0
	}

	// For the 'l' block:
	//
	//     W['d'] = sqrt( sk ./ zk )
	//     W['di'] = sqrt( zk ./ sk )
	//     lambdak = sqrt( sk .* zk )
	//
	// where sk and zk are the first dims['l'] entries of s and z.
	// lambda_k is stored in the first dims['l'] positions of lmbda.

	m := dims.At("l")[0]
	//td := s.FloatArray()
	stmp = matrix.FloatVector(s.FloatArray()[mnl : mnl+m])
	//zd := z.FloatArray()
	ztmp = matrix.FloatVector(z.FloatArray()[mnl : mnl+m])
	//fmt.Printf(".Sqrt()=\n%v\n", matrix.Div(stmp, ztmp).Sqrt().ToString("%.17f"))
	//d := stmp.Div(ztmp)
	//d.Apply(d, math.Sqrt)
	d := matrix.Div(stmp, ztmp).Sqrt()
	//di := d.Copy()
	//di.Apply(di, func(a float64)float64 { return 1.0/a })
	di := matrix.Inv(d)
	//fmt.Printf("d:\n%v\n", d)
	//fmt.Printf("di:\n%v\n", di)
	W.Set("d", d)
	W.Set("di", di)
	//lmd = stmp.Mul(ztmp)
	//lmd.Apply(lmd, math.Sqrt)
	lmd = matrix.Mul(stmp, ztmp).Sqrt()
	// lmd has indexes mnl:mnl+m and length of m
	lmbda.SetIndexesFromArray(lmd.FloatArray(), matrix.MakeIndexSet(mnl, mnl+m, 1)...)
	//fmt.Printf("after l:\n%v\n", lmbda)

	/*
	   For the 'q' blocks, compute lists 'v', 'beta'.

	   The vector v[k] has unit hyperbolic norm:

	       (sqrt( v[k]' * J * v[k] ) = 1 with J = [1, 0; 0, -I]).

	   beta[k] is a positive scalar.

	   The hyperbolic Householder matrix H = 2*v[k]*v[k]' - J
	   defined by v[k] satisfies

	       (beta[k] * H) * zk  = (beta[k] * H) \ sk = lambda_k

	   where sk = s[indq[k]:indq[k+1]], zk = z[indq[k]:indq[k+1]].

	   lambda_k is stored in lmbda[indq[k]:indq[k+1]].
	*/
	ind := mnl + dims.At("l")[0]
	var beta *matrix.FloatMatrix

	for _, k := range dims.At("q") {
		W.Append("v", matrix.FloatZeros(k, 1))
	}
	beta = matrix.FloatZeros(len(dims.At("q")), 1)
	W.Set("beta", beta)
	vset := W.At("v")
	for k, m := range dims.At("q") {
		v := vset[k]
		// a = sqrt( sk' * J * sk )  where J = [1, 0; 0, -I]
		aa := jnrm2(s, m, ind)
		// b = sqrt( zk' * J * zk )
		bb := jnrm2(z, m, ind)
		// beta[k] = ( a / b )**1/2
		beta.SetIndex(k, math.Sqrt(aa/bb))
		// c = sqrt( (sk/a)' * (zk/b) + 1 ) / sqrt(2)
		c0 := blas.DotFloat(s, z, &la_.IOpt{"n", m},
			&la_.IOpt{"offsetx", ind}, &la_.IOpt{"offsety", ind})
		cc := math.Sqrt((c0/aa/bb + 1.0) / 2.0)

		// vk = 1/(2*c) * ( (sk/a) + J * (zk/b) )
		blas.CopyFloat(z, v, &la_.IOpt{"offsetx", ind}, &la_.IOpt{"n", m})
		blas.ScalFloat(v, -1.0/bb)
		v.SetIndex(0, -1.0*v.GetIndex(0))
		blas.AxpyFloat(s, v, 1.0/aa, &la_.IOpt{"offsetx", ind}, &la_.IOpt{"n", m})
		blas.ScalFloat(v, 1.0/2.0/cc)

		// v[k] = 1/sqrt(2*(vk0 + 1)) * ( vk + e ),  e = [1; 0]
		v.SetIndex(0, v.GetIndex(0)+1.0)
		blas.ScalFloat(v, (1.0 / math.Sqrt(2.0*v.GetIndex(0))))
		/*
		   To get the scaled variable lambda_k

		       d =  sk0/a + zk0/b + 2*c
		       lambda_k = [ c;
		                    (c + zk0/b)/d * sk1/a + (c + sk0/a)/d * zk1/b ]
		       lambda_k *= sqrt(a * b)
		*/
		lmbda.SetIndex(ind, cc)
		dd := 2*cc + s.GetIndex(ind)/aa + z.GetIndex(ind)/bb
		blas.CopyFloat(s, lmbda, &la_.IOpt{"offsetx", ind + 1}, &la_.IOpt{"offsety", ind + 1},
			&la_.IOpt{"n", m - 1})
		zz := (cc + z.GetIndex(ind)/bb) / dd / aa
		ss := (cc + s.GetIndex(ind)/aa) / dd / bb
		blas.ScalFloat(lmbda, zz, &la_.IOpt{"offset", ind + 1}, &la_.IOpt{"n", m - 1})
		blas.AxpyFloat(z, lmbda, ss, &la_.IOpt{"offsetx", ind + 1},
			&la_.IOpt{"offsety", ind + 1}, &la_.IOpt{"n", m - 1})
		blas.ScalFloat(lmbda, math.Sqrt(aa*bb), &la_.IOpt{"offset", ind}, &la_.IOpt{"n", m})

		ind += m
		//fmt.Printf("after q[%d]:\n%v\n", k, lmbda)
	}
	/*
	   For the 's' blocks: compute two lists 'r' and 'rti'.

	       r[k]' * sk^{-1} * r[k] = diag(lambda_k)^{-1}
	       r[k]' * zk * r[k] = diag(lambda_k)

	   where sk and zk are the entries inds[k] : inds[k+1] of
	   s and z, reshaped into symmetric matrices.

	   rti[k] is the inverse of r[k]', so

	       rti[k]' * sk * rti[k] = diag(lambda_k)^{-1}
	       rti[k]' * zk^{-1} * rti[k] = diag(lambda_k).

	   The vectors lambda_k are stored in

	       lmbda[ dims['l'] + sum(dims['q']) : -1 ]
	*/
	for _, k := range dims.At("s") {
		W.Append("r", matrix.FloatZeros(k, k))
		W.Append("rti", matrix.FloatZeros(k, k))
	}
	maxs := maxdim(dims.At("s"))
	work := matrix.FloatZeros(maxs*maxs, 1)
	Ls := matrix.FloatZeros(maxs*maxs, 1)
	Lz := matrix.FloatZeros(maxs*maxs, 1)
	ind2 := ind
	for k, m := range dims.At("s") {
		r := W.At("r")[k]
		rti := W.At("rti")[k]

		// Factor sk = Ls*Ls'; store Ls in ds[inds[k]:inds[k+1]].
		blas.CopyFloat(s, Ls, &la_.IOpt{"offsetx", ind2}, &la_.IOpt{"n", m * m})
		lapack.PotrfFloat(Ls, &la_.IOpt{"n", m}, &la_.IOpt{"lda", m})

		// Factor zs[k] = Lz*Lz'; store Lz in dz[inds[k]:inds[k+1]].
		blas.CopyFloat(z, Lz, &la_.IOpt{"offsetx", ind2}, &la_.IOpt{"n", m * m})
		lapack.PotrfFloat(Lz, &la_.IOpt{"n", m}, &la_.IOpt{"lda", m})

		// SVD Lz'*Ls = U*diag(lambda_k)*V'.  Keep U in work.
		for i := 0; i < m; i++ {
			blas.ScalFloat(Ls, 0.0, &la_.IOpt{"offset", i * m}, &la_.IOpt{"n", i})
		}
		blas.CopyFloat(Ls, work, &la_.IOpt{"n", m * m})
		blas.TrmmFloat(Lz, work, 1.0, la_.OptTransA, &la_.IOpt{"lda", m}, &la_.IOpt{"ldb", m},
			&la_.IOpt{"n", m}, &la_.IOpt{"m", m})
		lapack.GesvdFloat(work, lmbda, nil, nil,
			la_.OptJobuO, &la_.IOpt{"lda", m}, &la_.IOpt{"offsetS", ind},
			&la_.IOpt{"n", m}, &la_.IOpt{"m", m})

		// r = Lz^{-T} * U
		blas.CopyFloat(work, r, &la_.IOpt{"n", m * m})
		blas.TrsmFloat(Lz, r, 1.0, la_.OptTransA,
			&la_.IOpt{"lda", m}, &la_.IOpt{"n", m}, &la_.IOpt{"m", m})

		// rti = Lz * U
		blas.CopyFloat(work, rti, &la_.IOpt{"n", m * m})
		blas.TrmmFloat(Lz, rti, 1.0,
			&la_.IOpt{"lda", m}, &la_.IOpt{"n", m}, &la_.IOpt{"m", m})

		// r := r * diag(sqrt(lambda_k))
		// rti := rti * diag(1 ./ sqrt(lambda_k))
		for i := 0; i < m; i++ {
			a := math.Sqrt(lmbda.GetIndex(ind + i))
			blas.ScalFloat(r, a, &la_.IOpt{"offset", m * i}, &la_.IOpt{"n", m})
			blas.ScalFloat(rti, 1.0/a, &la_.IOpt{"offset", m * i}, &la_.IOpt{"n", m})
		}
		ind += m
		ind2 += m * m
	}
	return
}

// Inner product of two vectors in S.
func sdot(x, y *matrix.FloatMatrix, dims *sets.DimensionSet, mnl int) float64 {
	/*DEBUGGED*/
	ind := mnl + dims.At("l")[0] + dims.Sum("q")
	a := blas.DotFloat(x, y, &la_.IOpt{"n", ind})
	for _, m := range dims.At("s") {
		a += blas.DotFloat(x, y, &la_.IOpt{"offsetx", ind}, &la_.IOpt{"offsety", ind},
			&la_.IOpt{"incx", m + 1}, &la_.IOpt{"incy", m + 1}, &la_.IOpt{"n", m})
		for j := 1; j < m; j++ {
			a += 2.0 * blas.DotFloat(x, y, &la_.IOpt{"offsetx", ind + j}, &la_.IOpt{"offsety", ind + j},
				&la_.IOpt{"incx", m + 1}, &la_.IOpt{"incy", m + 1}, &la_.IOpt{"n", m - j})
		}
		ind += m * m
	}
	return a
}

// Returns the norm of a vector in S
func snrm2(x *matrix.FloatMatrix, dims *sets.DimensionSet, mnl int) float64 {
	/*DEBUGGED*/
	return math.Sqrt(sdot(x, x, dims, mnl))
}

// Converts lower triangular matrix to symmetric.
// Fills in the upper triangular part of the symmetric matrix stored in
// x[offset : offset+n*n] using 'L' storage.
func symm(x *matrix.FloatMatrix, n, offset int) (err error) {
	/*DEBUGGED*/
	err = nil
	if n <= 1 {
		return
	}
	for i := 0; i < n; i++ {
		err = blas.Copy(x, x, &la_.IOpt{"offsetx", offset + i*(n+1) + 1},
			&la_.IOpt{"offsety", offset + (i+1)*(n+1) - 1}, &la_.IOpt{"incy", n},
			&la_.IOpt{"n", n - i - 1})
		//if err != nil { return }
	}
	return
}

/*
   Matrix-vector multiplication.

   A is a matrix or spmatrix of size (m, n) where

       N = dims['l'] + sum(dims['q']) + sum( k**2 for k in dims['s'] )

   representing a mapping from R^n to S.

   If trans is 'N':

       y := alpha*A*x + beta * y   (trans = 'N').

   x is a vector of length n.  y is a vector of length N.

   If trans is 'T':

       y := alpha*A'*x + beta * y  (trans = 'T').

   x is a vector of length N.  y is a vector of length n.

   The 's' components in S are stored in unpacked 'L' storage.
*/
func sgemv(A, x, y *matrix.FloatMatrix, alpha, beta float64, dims *sets.DimensionSet, opts ...la_.Option) error {

	m := dims.Sum("l", "q") + dims.SumSquared("s")
	n := la_.GetIntOpt("n", -1, opts...)
	if n == -1 {
		n = A.Cols()
	}
	trans := la_.GetIntOpt("trans", int(la_.PNoTrans), opts...)
	offsetX := la_.GetIntOpt("offsetx", 0, opts...)
	offsetY := la_.GetIntOpt("offsety", 0, opts...)
	offsetA := la_.GetIntOpt("offseta", 0, opts...)

	if trans == int(la_.PTrans) && alpha != 0.0 {
		trisc(x, dims, offsetX)
		//fmt.Printf("trisc x=\n%v\n", x.ConvertToString())
	}
	//fmt.Printf("alpha=%.4f beta=%.4f m=%d n=%d\n", alpha, beta, m, n)
	//fmt.Printf("A=\n%v\nx=\n%v\ny=\n%v\n", A, x.ConvertToString(), y.ConvertToString())
	err := blas.GemvFloat(A, x, y, alpha, beta, &la_.IOpt{"trans", trans},
		&la_.IOpt{"n", n}, &la_.IOpt{"m", m}, &la_.IOpt{"offseta", offsetA},
		&la_.IOpt{"offsetx", offsetX}, &la_.IOpt{"offsety", offsetY})
	//fmt.Printf("gemv y=\n%v\n", y.ConvertToString())

	if trans == int(la_.PTrans) && alpha != 0.0 {
		triusc(x, dims, offsetX)
	}
	return err
}

/*
 The inverse product x := (y o\ x), when the 's' components of y are
 diagonal.
*/

func sinv(x, y *matrix.FloatMatrix, dims *sets.DimensionSet, mnl int) (err error) {
	/*DEBUGGED*/

	err = nil

	// For the nonlinear and 'l' blocks:
	//
	//     yk o\ xk = yk .\ xk.

	ind := mnl + dims.At("l")[0]
	blas.Tbsv(y, x, &la_.IOpt{"n", ind}, &la_.IOpt{"k", 0}, &la_.IOpt{"ldA", 1})

	// For the 'q' blocks:
	//
	//                        [ l0   -l1'              ]
	//     yk o\ xk = 1/a^2 * [                        ] * xk
	//                        [ -l1  (a*I + l1*l1')/l0 ]
	//
	// where yk = (l0, l1) and a = l0^2 - l1'*l1.

	for _, m := range dims.At("q") {
		aa := blas.Nrm2Float(y, &la_.IOpt{"n", m - 1}, &la_.IOpt{"offset", ind + 1})
		ee := y.GetIndex(ind)
		aa = (ee + aa) * (ee - aa)
		cc := x.GetIndex(ind)
		dd := blas.DotFloat(x, y, &la_.IOpt{"n", m - 1}, &la_.IOpt{"offsetx", ind + 1},
			&la_.IOpt{"offsety", ind + 1})
		x.SetIndex(ind, cc*ee-dd)
		blas.ScalFloat(x, aa/ee, &la_.IOpt{"n", m - 1}, &la_.IOpt{"offset", ind + 1})
		blas.AxpyFloat(y, x, dd/ee-cc, &la_.IOpt{"n", m - 1},
			&la_.IOpt{"offsetx", ind + 1}, &la_.IOpt{"offsety", ind + 1})
		blas.ScalFloat(x, 1.0/aa, &la_.IOpt{"n", m}, &la_.IOpt{"offset", ind})
		ind += m
	}

	// For the 's' blocks:
	//
	//     yk o\ xk =  xk ./ gamma
	//
	// where gammaij = .5 * (yk_i + yk_j).

	ind2 := ind
	for _, m := range dims.At("s") {
		for j := 0; j < m; j++ {
			u := matrix.FloatVector(y.FloatArray()[ind2+j : ind2+m])
			u.Add(y.GetIndex(ind2 + j))
			u.Scale(0.5)
			blas.Tbsv(u, x, &la_.IOpt{"n", m - j}, &la_.IOpt{"k", 0}, &la_.IOpt{"lda", 1},
				&la_.IOpt{"offsetx", ind + j*(m+1)})
		}
		ind += m * m
		ind2 += m
	}
	return
}

/*
   Sets upper triangular part of the 's' components of x equal to zero
   and scales the strictly lower triangular part by 2.0.
*/
func trisc(x *matrix.FloatMatrix, dims *sets.DimensionSet, offset int) error {

	//m := dims.Sum("l", "q") + dims.SumSquared("s")
	ind := offset + dims.Sum("l", "q")

	for _, mk := range dims.At("s") {
		for j := 1; j < mk; j++ {
			blas.ScalFloat(x, 0.0, la_.IntOpt("n", mk-j), la_.IntOpt("inc", mk),
				la_.IntOpt("offset", ind+j*(mk+1)-1))
			blas.ScalFloat(x, 2.0, la_.IntOpt("n", mk-j), la_.IntOpt("offset", ind+mk*(j-1)+j))
		}
		ind += mk * mk
	}
	return nil
}

/*
   Scales the strictly lower triangular part of the 's' components of x
   by 0.5.

*/
func triusc(x *matrix.FloatMatrix, dims *sets.DimensionSet, offset int) error {

	//m := dims.Sum("l", "q") + dims.SumSquared("s")
	ind := offset + dims.Sum("l", "q")

	for _, mk := range dims.At("s") {
		for j := 1; j < mk; j++ {
			blas.ScalFloat(x, 0.5, &la_.IOpt{"n", mk - j}, &la_.IOpt{"offset", ind + mk*(j-1) + j})
		}
		ind += mk * mk
	}
	return nil
}

// The product x := (y o x).  If diag is 'D', the 's' part of y is
// diagonal and only the diagonal is stored.
func sprod(x, y *matrix.FloatMatrix, dims *sets.DimensionSet, mnl int, opts ...la_.Option) (err error) {

	err = nil
	diag := la_.GetStringOpt("diag", "N", opts...)
	// For the nonlinear and 'l' blocks:
	//
	//     yk o xk = yk .* xk.
	ind := mnl + dims.At("l")[0]
	err = blas.Tbmv(y, x, &la_.IOpt{"n", ind}, &la_.IOpt{"k", 0}, &la_.IOpt{"lda", 1})
	if err != nil {
		return
	}
	//fmt.Printf("Sprod l:x=\n%v\n", x)

	// For 'q' blocks:
	//
	//               [ l0   l1'  ]
	//     yk o xk = [           ] * xk
	//               [ l1   l0*I ]
	//
	// where yk = (l0, l1).
	for _, m := range dims.At("q") {
		dd := blas.DotFloat(x, y, &la_.IOpt{"offsetx", ind}, &la_.IOpt{"offsety", ind},
			&la_.IOpt{"n", m})
		//fmt.Printf("dd=%v\n", dd)
		alpha := y.GetIndex(ind)
		//fmt.Printf("scal=%v\n", alpha)
		blas.ScalFloat(x, alpha, &la_.IOpt{"offset", ind + 1}, &la_.IOpt{"n", m - 1})
		alpha = x.GetIndex(ind)
		//fmt.Printf("axpy=%v\n", alpha)
		blas.AxpyFloat(y, x, alpha, &la_.IOpt{"offsetx", ind + 1}, &la_.IOpt{"offsety", ind + 1},
			&la_.IOpt{"n", m - 1})
		x.SetIndex(ind, dd)
		ind += m
	}
	//fmt.Printf("Sprod q :x=\n%v\n", x)

	// For the 's' blocks:
	//
	//    yk o sk = .5 * ( Yk * mat(xk) + mat(xk) * Yk )
	//
	// where Yk = mat(yk) if diag is 'N' and Yk = diag(yk) if diag is 'D'.

	if diag[0] == 'N' {
		// DEBUGGED
		maxm := maxdim(dims.At("s"))
		A := matrix.FloatZeros(maxm, maxm)
		for _, m := range dims.At("s") {
			blas.Copy(x, A, &la_.IOpt{"offsetx", ind}, &la_.IOpt{"n", m * m})
			for i := 0; i < m-1; i++ { // i < m-1 --> i < m
				symm(A, m, 0)
				symm(y, m, ind)
			}
			err = blas.Syr2kFloat(A, y, x, 0.5, 0.0, &la_.IOpt{"n", m}, &la_.IOpt{"k", m},
				&la_.IOpt{"lda", m}, &la_.IOpt{"ldb", m}, &la_.IOpt{"ldc", m},
				&la_.IOpt{"offsetb", ind}, &la_.IOpt{"offsetc", ind})
			if err != nil {
				return
			}
			ind += m * m
		}
		//fmt.Printf("Sprod diag=N s:x=\n%v\n", x)

	} else {
		ind2 := ind
		for _, m := range dims.At("s") {
			for i := 0; i < m; i++ {
				// original: u = 0.5 * ( y[ind2+i:ind2+m] + y[ind2+i] )
				// creates matrix of elements: [ind2+i ... ind2+m] then
				// element wisely adds y[ind2+i] and scales by 0.5
				iset := matrix.MakeIndexSet(ind2+i, ind2+m, 1)
				u := matrix.FloatVector(y.GetIndexes(iset...))
				u.Add(y.GetIndex(ind2 + i))
				u.Scale(0.5)
				err = blas.Tbmv(u, x, &la_.IOpt{"n", m - i}, &la_.IOpt{"k", 0}, &la_.IOpt{"lda", 1},
					&la_.IOpt{"offsetx", ind + i*(m+1)})
				if err != nil {
					return
				}
			}
			ind += m * m
			ind2 += m
		}
		//fmt.Printf("Sprod diag=T s:x=\n%v\n", x)
	}
	return
}

// The product x := y o y.   The 's' components of y are diagonal and
// only the diagonals of x and y are stored.
func ssqr(x, y *matrix.FloatMatrix, dims *sets.DimensionSet, mnl int) (err error) {
	/*DEBUGGED*/
	blas.Copy(y, x)
	ind := mnl + dims.At("l")[0]
	err = blas.Tbmv(y, x, &la_.IOpt{"n", ind}, &la_.IOpt{"k", 0}, &la_.IOpt{"lda", 1})
	if err != nil {
		return
	}

	for _, m := range dims.At("q") {
		v := blas.Nrm2Float(y, &la_.IOpt{"n", m}, &la_.IOpt{"offset", ind})
		x.SetIndex(ind, v*v)
		blas.ScalFloat(x, 2.0*y.GetIndex(ind), &la_.IOpt{"n", m - 1}, &la_.IOpt{"offset", ind + 1})
		ind += m
	}
	err = blas.Tbmv(y, x, &la_.IOpt{"n", dims.Sum("s")}, &la_.IOpt{"k", 0},
		&la_.IOpt{"lda", 1}, &la_.IOpt{"offseta", ind}, &la_.IOpt{"offsetx", ind})
	return
}

// Returns min {t | x + t*e >= 0}, where e is defined as follows
//
//  - For the nonlinear and 'l' blocks: e is the vector of ones.
//  - For the 'q' blocks: e is the first unit vector.
//  - For the 's' blocks: e is the identity matrix.
//
// When called with the argument sigma, also returns the eigenvalues
// (in sigma) and the eigenvectors (in x) of the 's' components of x.
func maxStep(x *matrix.FloatMatrix, dims *sets.DimensionSet, mnl int, sigma *matrix.FloatMatrix) (rval float64, err error) {
	/*DEBUGGED*/

	rval = 0.0
	err = nil
	t := make([]float64, 0, 10)
	ind := mnl + dims.Sum("l")
	if ind > 0 {
		t = append(t, -minvec(x.FloatArray()[:ind]))
	}
	for _, m := range dims.At("q") {
		if m > 0 {
			v := blas.Nrm2Float(x, &la_.IOpt{"offset", ind + 1}, &la_.IOpt{"n", m - 1})
			v -= x.GetIndex(ind)
			t = append(t, v)
		}
		ind += m
	}

	//var Q *matrix.FloatMatrix
	//var w *matrix.FloatMatrix
	ind2 := 0
	//if sigma == nil && len(dims.At("s")) > 0 {
	//	mx := dims.Max("s")
	//	Q = matrix.FloatZeros(mx, mx)
	//	w = matrix.FloatZeros(mx, 1)
	//}
	for _, m := range dims.At("s") {
		if sigma == nil {
			Q := matrix.FloatZeros(m, m)
			w := matrix.FloatZeros(m, 1)
			blas.Copy(x, Q, &la_.IOpt{"offsetx", ind}, &la_.IOpt{"n", m * m})
			err = lapack.SyevrFloat(Q, w, nil, 0.0, nil, []int{1, 1}, la_.OptRangeInt,
				&la_.IOpt{"n", m}, &la_.IOpt{"lda", m})
			if m > 0 && err == nil {
				t = append(t, -w.GetIndex(0))
			}
		} else {
			err = lapack.SyevdFloat(x, sigma, la_.OptJobZValue, &la_.IOpt{"n", m},
				&la_.IOpt{"lda", m}, &la_.IOpt{"offseta", ind}, &la_.IOpt{"offsetw", ind2})
			if m > 0 {
				t = append(t, -sigma.GetIndex(ind2))
			}
		}
		ind += m * m
		ind2 += m
	}

	if len(t) > 0 {
		rval = maxvec(t)
	}
	return
}

/*
   Copy x to y using packed storage.

   The vector x is an element of S, with the 's' components stored in
   unpacked storage.  On return, x is copied to y with the 's' components
   stored in packed storage and the off-diagonal entries scaled by
   sqrt(2).
*/
func pack(x, y *matrix.FloatMatrix, dims *sets.DimensionSet, opts ...la_.Option) (err error) {
	/*DEBUGGED*/
	err = nil
	mnl := la_.GetIntOpt("mnl", 0, opts...)
	offsetx := la_.GetIntOpt("offsetx", 0, opts...)
	offsety := la_.GetIntOpt("offsety", 0, opts...)

	nlq := mnl + dims.At("l")[0] + dims.Sum("q")
	blas.Copy(x, y, &la_.IOpt{"n", nlq}, &la_.IOpt{"offsetx", offsetx},
		&la_.IOpt{"offsety", offsety})

	iu, ip := offsetx+nlq, offsety+nlq
	for _, n := range dims.At("s") {
		for k := 0; k < n; k++ {
			blas.Copy(x, y, &la_.IOpt{"n", n - k}, &la_.IOpt{"offsetx", iu + k*(n+1)},
				&la_.IOpt{"offsety", ip})
			y.SetIndex(ip, (y.GetIndex(ip) / math.Sqrt(2.0)))
			ip += n - k
		}
		iu += n * n
	}
	np := dims.SumPacked("s")
	blas.ScalFloat(y, math.Sqrt(2.0), &la_.IOpt{"n", np}, &la_.IOpt{"offset", offsety + nlq})
	return
}

// In-place version of pack(), which also accepts matrix arguments x.
// The columns of x are elements of S, with the 's' components stored
// in unpacked storage.  On return, the 's' components are stored in
// packed storage and the off-diagonal entries are scaled by sqrt(2).
//
func pack2(x *matrix.FloatMatrix, dims *sets.DimensionSet, mnl int) (err error) {
	if len(dims.At("s")) == 0 {
		return nil
	}

	const sqrt2 = 1.41421356237309504880

	iu := mnl + dims.Sum("l", "q")
	ip := iu
	row := matrix.FloatZeros(1, x.Cols())
	//fmt.Printf("x.size = %d %d\n", x.Rows(), x.Cols())
	for _, n := range dims.At("s") {
		for k := 0; k < n; k++ {
			cnt := n - k
			row = x.GetRow(iu+(n+1)*k, row)
			//fmt.Printf("%02d: %v\n", iu+(n+1)*k, x.FloatArray())
			x.SetRow(ip, row)
			for i := 1; i < n-k; i++ {
				row = x.GetRow(iu+(n+1)*k+i, row)
				//fmt.Printf("%02d: %v\n", iu+(n+1)*k+i, x.FloatArray())
				x.SetRow(ip+i, row.Scale(sqrt2))
			}
			ip += cnt
		}
		iu += n * n
	}
	return nil
}

/*
   The vector x is an element of S, with the 's' components stored
   in unpacked storage and off-diagonal entries scaled by sqrt(2).
   On return, x is copied to y with the 's' components stored in
   unpacked storage.

*/
func unpack(x, y *matrix.FloatMatrix, dims *sets.DimensionSet, opts ...la_.Option) (err error) {
	/*DEBUGGED*/
	err = nil
	mnl := la_.GetIntOpt("mnl", 0, opts...)
	offsetx := la_.GetIntOpt("offsetx", 0, opts...)
	offsety := la_.GetIntOpt("offsety", 0, opts...)

	nlq := mnl + dims.At("l")[0] + dims.Sum("q")
	err = blas.Copy(x, y, &la_.IOpt{"n", nlq}, &la_.IOpt{"offsetx", offsetx},
		&la_.IOpt{"offsety", offsety})
	if err != nil {
		return
	}

	ip, iu := offsetx+nlq, offsety+nlq
	for _, n := range dims.At("s") {
		for k := 0; k < n; k++ {
			err = blas.Copy(x, y, &la_.IOpt{"n", n - k}, &la_.IOpt{"offsetx", ip},
				&la_.IOpt{"offsety", iu + k*(n+1)})
			if err != nil {
				return
			}

			ip += n - k
			blas.ScalFloat(y, 1.0/math.Sqrt(2.0),
				&la_.IOpt{"n", n - k - 1}, &la_.IOpt{"offset", iu + k*(n+1) + 1})
		}
		iu += n * n
	}
	/*
		nu := dims.SumSquared("s")
		fmt.Printf("-- UnPack: nu=%d, offset=%d\n", nu, offsety+nlq)
		err = blas.ScalFloat(y,
			&la_.IOpt{"n", nu}, &la_.IOpt{"offset", offsety+nlq})
	*/
	return
}

/*
   Returns x' * J * y, where J = [1, 0; 0, -I].
*/
func jdot(x, y *matrix.FloatMatrix, n, offsetx, offsety int) float64 {
	if n <= 0 {
		n = x.NumElements()
	}
	a := blas.DotFloat(x, y, &la_.IOpt{"n", n - 1}, &la_.IOpt{"offsetx", offsetx + 1},
		&la_.IOpt{"offsety", offsety + 1})
	return x.GetIndex(offsetx)*y.GetIndex(offsety) - a
}

/*
   Returns sqrt(x' * J * x) where J = [1, 0; 0, -I], for a vector
   x in a second order cone.
*/
func jnrm2(x *matrix.FloatMatrix, n, offset int) float64 {
	/*DEBUGGED*/
	if n <= 0 {
		n = x.NumElements()
	}
	if offset < 0 {
		offset = 0
	}
	a := blas.Nrm2Float(x, &la_.IOpt{"n", n - 1}, &la_.IOpt{"offset", offset + 1})
	fst := x.GetIndex(offset)
	return math.Sqrt(fst-a) * math.Sqrt(fst+a)
}

// Local Variables:
// tab-width: 4
// End:
