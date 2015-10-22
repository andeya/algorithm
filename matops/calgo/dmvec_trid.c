
// Copyright (c) Harri Rautila, 2012,2013

// This file is part of github.com/hrautila/matops package. It is free software,
// distributed under the terms of GNU Lesser General Public License Version 3, or
// any later version. See the COPYING tile included in this archive.

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "cmops.h"
#include "inner_vec_axpy.h"
#include "inner_vec_dot.h"

/*
 TRMV UPPER:

     A00 : a01 | A02       x0
     ---------------      ----   
      0  : a11 | a12       x1
     ===============      ====
      0  :  0  | A22       x2

   x0 = A00*x0 + a01*x1 + A02*x2
   x1 =          a11*x1 + a12*x2
   x2 =                   A22*x2

   x0 not needed for x2; start from top breadth first; update x0 with current
   column [a01; a11], then update x1

   AXPY version: (axpy_forward)
     x0 = x0 + a01*x1
     x1 = a11*x1 

   UPPER, TRANSA

     A00 |  0  :  0        x0
     ===============      ====   
     a10 | a11 :  0        x1
     ---------------      ----
     A20 | a12 : A22       x2

   A00 = A00, a01 == a10, A02 == A20, a12 == a21
   
   x0 not need for x2; do depth first from bottom to top.

   DOT version: (ddot_backward)
     x1 = a11*x11 + a12*x2.T


 TRMV LOWER

     A00 |  0  :  0        x0
     ===============      ====        Xt
     a10 | a11 :  0        x1   ---> ----
     ---------------      ----        Xl
     A20 | a12 : A22       x2

     x0 = A00*x0
     x1 = a10*x0 + a11*x11 
     x2 = A20*x0 + a12*x11 + A22*x2

     x2 not needed for x1, start backward from bottom
     
     AXPY version (axpy_backward)
       x2 = x2 + a12*x2.T
       x1 = a11*x11

 */

// Functions here implement various versions of TRMV operation.


/*

 TRMV LOWER, NOTRANS

     A00 |  0  :  0        x0
     ===============      ====  
     a10 | a11 :  0        x1   
     ---------------      ----  
     A20 | a12 : A22       x2

   x0 = A00*x0
   x1 = a10*x0 + a11*x1
   x2 = A20*x0 + a21*x1 + A22*x2

   x2 not needed for x1; start backward from bottom

   Calculates backward a diagonal block and updates Xc values from last to first.
   Updates are calculated in breadth first manner by with successive AXPY operations.
 */
static void
_dmvec_trid_axpy_backward(double *Xc, const double *Ac, int unit,
                          int incX, int ldA, int nRE)
{
  // Y is 
  register int i;
  register double *xr, xtmp;
  register const double *Ar, *Acl;

  // diagonal matrix of nRE rows/cols and vector X of length nRE
  // move to point to last column and last entry of X.
  Acl = Ac + (nRE-1)*ldA;
  xr = Xc + (nRE-1)*incX;

  // xr is the current X element, Ar is row in current A column.
  for (i = nRE; i > 0; i--) {
    Ar = Acl + i - 1; // move on diagonal

    // update all x-values below with the current A column and current X
    _inner_vec_daxpy(xr+incX, incX, Ar+1, xr, incX, 1.0, nRE-i);
    xr[0] *= unit ? 1.0 : Ar[0];

    // previous X, previous column in A 
    xr  -= incX;
    Acl -= ldA;
  }
}

// Calculates backward a diagonal block and updates Xc values from bottom to top.
// Updates are calculated in depth first manner with DOT operations.
static void
_dmvec_trid_dot_backward(double *Xc, const double *Ac, int unit,
                         int incX, int ldA, int nRE)
{
  // Y is 
  register int i;
  register double *xr;
  double xtmp;
  register const double *Ar, *Acl;

  // lower diagonal matrix (transposed) of nRE rows/cols and vector X of length nRE
  // we do it really forward!! unlike the _axpy method above.
  Acl = Ac + (nRE-1)*ldA;
  xr = Xc + (nRE-1)*incX;

  // xr is the current X element, Ar is row in current A column.
  for (i = 0; i < nRE; i++) {
    //Ar = Ac + i; // move on diagonal

    // update current x-value with the current A column and top
    xtmp = unit ? xr[0] : 0.0;
    _inner_vec_ddot(&xtmp, 1, Acl, Xc, incX, 1.0, nRE-unit-i);
    xr[0] = xtmp;

    // previous X, previous column in A 
    xr -= incX;
    Acl -= ldA;
  }
}

// Calculate forward a diagonal block and updates x0 values from first to last.
// x0 updated in breadth first manner with successive AXPY operations.
static void
_dmvec_trid_axpy_forward(double *Xc, const double *Ac, double unit,
                         int incX, int ldA, int nRE)
{
  // Y is 
  register int i;
  register double *X0, *x1;
  register const double *a11;

  // upper diagonal matrix of nRE rows/cols and vector X, Y of length nRE
  X0 = Xc;
  x1 = Xc;

  // xr is the current X element, Ar is row in current A column.
  for (i = 0; i < nRE; i++) {
    // update all previous x-values with current A column and current X
    _inner_vec_daxpy(X0, incX, Ac, x1, incX, 1.0, i);
    a11 = Ac + i;
    x1[0] *= unit ? 1.0 : a11[0];
    // next X, next column in A 
    x1 += incX;
    Ac += ldA;
  }
}

// Update X values from top to bottom; current x is inner product of current A [a11; a12]
// column and current and following x values [x1; x2]
static void
_dmvec_trid_dot_forward(double *Xc, const double *Ac, int unit, int incX, int ldA, int nRE)
{
  // Y is 
  register int i;
  register double *x1;
  double xtmp;
  register const double *a11;

  x1 = Xc;
  // xr is the current X element, Ar is row in current A column.
  for (i = 0; i < nRE; i++) {
    a11 = Ac + i + unit;
    // update current x-value with current A column and current and following X
    xtmp = unit ? x1[0] : 0.0;
    _inner_vec_ddot(&xtmp, 1, a11, x1, incX, 1.0, nRE-unit-i);
    x1[0] = xtmp;
    // next X, next column in A 
    x1 += incX;
    Ac += ldA;
  }
}

//extern void memset(void *, int, size_t);

// X = A*X; unblocked version
void dmvec_trid_unb(mvec_t *X, const mdata_t *A, int flags, int N)
{
  // indicates if diagonal entry is unit (=1.0) or non-unit.
  int unit = flags & MTX_UNIT ? 1 : 0;

  if (N <= 0) {
    return;
  }
  if (flags & MTX_UPPER) {
    if (flags & MTX_TRANSA || flags & MTX_TRANS) {
      _dmvec_trid_dot_backward(X->md, A->md, unit, X->inc, A->step, N);
    } else {
      _dmvec_trid_axpy_forward(X->md, A->md, unit, X->inc, A->step, N);
    }
  } else {
    if (flags & MTX_TRANSA || flags & MTX_TRANS) {
      _dmvec_trid_dot_forward(X->md, A->md, unit, X->inc, A->step, N);
    } else {
      _dmvec_trid_axpy_backward(X->md, A->md, unit, X->inc, A->step, N);
    }
  }
}


void dmvec_trid_blocked(mvec_t *X, const mdata_t *A, double alpha, int flags, int N, int NB)
{
  int i, nI;
  mvec_t Y;
  // indicates if diagonal entry is unit (=1.0) or non-unit.
  int unit = flags & MTX_UNIT ? 1 : 0;

  if (N <= 0) {
    return;
  }
  if (NB <= 0) {
    NB = 68;
  }
  //memset(cB, 0, sizeof(cB));

  if (flags & MTX_UPPER) {
    for (i = 0; i < N; i += NB) {
      nI = N - i < NB ? N - i : NB;
      _dmvec_trid_axpy_forward(&X->md[i], &A->md[i*A->step+i], unit, X->inc, A->step, nI);
    }
  }
}

// Local Variables:
// indent-tabs-mode: nil
// End:
