
// Copyright (c) Harri Rautila, 2012,2013

// This file is part of github.com/hrautila/matops package. It is free software,
// distributed under the terms of GNU Lesser General Public License Version 3, or
// any later version. See the COPYING tile included in this archive.

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "cmops.h"
#include "inner_vec_axpy.h"

/*
  A: N*N, lower          X: N*1
     A00 | a01 | A02       x0
     ---------------      ----   
     a10 | a11 | a12       x1
     ---------------      ----
     A20 | a21 | A22       x2

   i = n: dimensions,
       a11 = 1*1, a10 = 1*n, a01 = n*1, a12 = N-(n+1)*1, A00 = n*n
       x1  = 1*1, x0  = n*1, x2  = N-(n+1)*1  

   A02, a12 are zeros, A00, A22 are lower tridiagonal

   Variant 0:
      if A diagonal is non-UNIT
         x1 = (x1 - a10*x0) / a11
      else
         x1 = x1 - a10*x0
   
   Variant 1:
      if A diagonal is non-UNIT
         x1 = x1/a11
         x2 = x2 - a21*x1;
      else
         x1 = x1
         x2 = x2 - a21*x1

   Variant 0 involves operation a10*x0 that is a DOT operation, and results 2n memory reads
   and 1 write. Vector a10 is row vector with elements N elements a part and x0 is column vector.
   Accessing x[i+1] is moving in memory direction and is likely a cache-hit. Accesing
   a10[i+1] is a memory reference to N*sizeof(double) bytes far and a likely cache-miss.
   
   Variant 1 update x2 with AXPY operation. Vector a21 is column vector and access to
   a21[i+1] is likely a cache-hit. A drawback is that we have in addition to read x2
   also write x2 which all in all causes N*N/2 read and write operations to x.

 */

// solves backward a diagonal block and updates Xc values.  Yc vector
// contains per row accumulated values for already solved X's. 
static void
_dmvec_solve_backward(double *Xc, const double *Ac, int unit, 
                      int incX, int ldA, int nRE)
{
  // Y is 
  register int i;
  double *x0, *x1;
  const double *a11, *A01;

  // upper diagonal matrix of nRE rows/cols and vector X, Y of length nRE
  // move to point to last column and last entry of X.
  A01 = Ac + (nRE-1)*ldA;
  x1 = Xc + (nRE-1)*incX;
  x0 = Xc;

  // xr is the current X element, Ar is row in current A column.
  for (i = nRE-1; i >= 0; i--) {
    a11 = A01 + i;                // move on diagonal
    x1[0] = unit ? x1[0] : x1[0]/a11[0];

    // update all x0-values with in current column (i is the count above current row)
    _inner_vec_daxpy(x0, incX, A01, x1, incX, -1.0, i);
    // repartition: previous X, previous column in A 
    x1  -= incX;
    A01 -= ldA;
  }
}

// solves forward a diagonal block and updates Xc values.  Yc vector
// contains per row accumulated values for already solved X's. 
static void
_dmvec_solve_forward(double *Xc, const double *Ac, int unit,
                     int incX, int ldA, int nRE)
{
  // Y is 
  register int i;
  register double *x1, *x2, xtmp;
  const double *a11, *a21;

  // lower diagonal matrix of nRE rows/cols and vector X, Y of length nRE
  a11 = Ac;
  x1 = Xc;
  x2 = x1 + incX;

  // xr is the current X element, Ar is row in current A column.
  for (i = 0; i < nRE; i++) {
    a11 = Ac + i;                // move on diagonal
    a21 = a11 + 1;
    //xtmp = xr[0] - yr[0];
    x1[0] = unit ? x1[0] : x1[0]/a11[0];
    // update all x2-values with in current column
    _inner_vec_daxpy(x2, incX, a21, x1, incX, -1.0, nRE-1-i);
    // next X, next column in A 
    x1 += incX;
    x2 = x1 + incX;
    Ac += ldA;
  }
}

extern void memset(void *, int, size_t);

#define MAX_VEC_NB 256

// X = A(-1)*X
void dmvec_solve_unb(mvec_t *X, const mdata_t *A, int flags, int N)
{
  int i, nI;
  int unit = flags & MTX_UNIT ? 1 : 0;

  if (flags & MTX_LOWER) {
    _dmvec_solve_forward(X->md, A->md, unit, X->inc, A->step, N);
  } else {
    _dmvec_solve_backward(X->md, A->md, unit, X->inc, A->step, N);
  }
}

/*
  A: N*N, lower          X: N*1

     A00 :  0  |  0        X0 
     ---------------      ----
     A10 : A11 |  0        X1 
     ===============      ====
     A20 : A21 | A22       X2 

   i = k: dimensions,
       A11 = n*n, A10 = n*k, A01 = k*n, A12 = N-(n+k)*n, A00 = k*k
       X1  = n*1, X0  = k*1, X2  = N-(n+k)*1  

   Lower Variant 0:
        X1 = X1 - A10*X0
        solve_forward_unb(X1, A11)

   Lower Variant 2:
        solve_forward_unb(X1, A11)
        X2 = X2 - A21 * X1

   UPPER

     A00 | A01 : A02       X0 
     ===============      ====
      0  | A11 : A12       X1 
     ---------------      ----
      0  |  0  : A22       X2 

   Lower Variant 2:
        solve_forward_unb(X1, A11)
        X0 = X0 - A01 * X1
*/
void dmvec_solve_blocked(mvec_t *X, const mdata_t *A, int flags, int N, int NB)
{
  int i, nI, nR;
  mdata_t A01, A11, A21;
  mvec_t X0, X1, X2;
  int unit = flags & MTX_UNIT ? 1 : 0;

  if (NB <= 0) {
    NB = 68;
  }

  A01.step = A->step;
  A11.step = A->step;
  A21.step = A->step;
  X0.inc = X->inc;
  X1.inc = X->inc;
  X2.inc = X->inc;
  X0.md = X->md;

  if (flags & MTX_LOWER) {
    nR = 0;
    for (i = 0; i < N; i += NB) {
      nI = N - i < NB ? N - i : NB;

      // Solve forward curent block
      X1.md = &X->md[i*X->inc];
      A11.md = &A->md[i*A->step+i];
      _dmvec_solve_forward(X1.md, A11.md, unit, X->inc, A->step, nI);
      nR += nI;

      // update X2 with new solutions.
      A21.md = &A->md[i*A->step+nR];
      X2.md = &X->md[nR*X->inc];
      //printf("A21:\n"); print_tile(A21.md, A21.step, N-nR, nI);
      dmult_gemv_blocked(&X2, &A21, &X1, -1.0, 1.0, 0, 0, nI, 0, N-nR, 0, 0);
    }
  } else {
    int n;
    nR = N;
    for (i = N; i > 0; i -= NB) {
      nI = i < NB ? i : NB;
      n  = i < NB? 0 : i-NB;

      // Solve forward curent block
      X1.md = &X->md[n*X->inc];
      A11.md = &A->md[n*A->step+n];
      _dmvec_solve_backward(X1.md, A11.md, unit, X->inc, A->step, nI);
      nR -= nI;

      // update X0 with new solutions.
      A01.md = &A->md[n*A->step];
      //printf("A01:\n"); print_tile(A01.md, A21.step, nR, nI);
      dmult_gemv_blocked(&X0, &A01, &X1, -1.0, 1.0, 0, 0, nI, 0, nR, 0, 0);
    }
  }
}

// Local Variables:
// indent-tabs-mode: nil
// End:
