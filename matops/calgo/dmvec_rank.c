
// Copyright (c) Harri Rautila, 2012,2013

// This file is part of github.com/hrautila/matops package. It is free software,
// distributed under the terms of GNU Lesser General Public License Version 3, or
// any later version. See the COPYING tile included in this archive.

#include <stdio.h>
#include <stdint.h>

#include "cmops.h"

// update 1 column of A matrix (a0) with vector X scaled with elements y0
static inline
void _inner_mv_daxpy(double *a0, const double *x, int incX,
		     const double *y0, double alpha, int nRE)
{
  register int i;
  register double cf0, cf1;

  cf0 = alpha * y0[0];
  
  for (i = 0; i < nRE-3; i += 4) {
    a0[0] += x[0] * cf0;
    x += incX;
    a0[1] += x[0] * cf0;
    x += incX;
    a0[2] += x[0] * cf0;
    x += incX;
    a0[3] += x[0] * cf0;
    x += incX;
    a0 += 4;
  }
  if (i == nRE)
    return;

  if (i < nRE-1) {
    a0[0] += x[0] * cf0;
    x += incX;
    a0[1] += x[0] * cf0;
    x += incX;
    a0 += 2;
    i += 2;
  }
  if (i < nRE) {
    a0[0] += x[0] * cf0;
    x += incX;
  }
}

// update 2 columns of A matrix (a0, a1) with vector X scaled with elements y0, y1
static inline
void _inner_mv_daxpy2(double *a0, double *a1, const double *x, int incX,
		      const double *y0, const double *y1, double alpha, int nRE)
{
  register int i;
  register double cf0, cf1;

  cf0 = alpha * y0[0];
  cf1 = alpha * y1[0];
  
  for (i = 0; i < nRE-3; i += 4) {
    a0[0] += x[0] * cf0;
    a1[0] += x[0] * cf1;
    x += incX;
    a0[1] += x[0] * cf0;
    a1[1] += x[0] * cf1;
    x += incX;
    a0[2] += x[0] * cf0;
    a1[2] += x[0] * cf1;
    x += incX;
    a0[3] += x[0] * cf0;
    a1[3] += x[0] * cf1;
    x += incX;
    a0 += 4;
    a1 += 4;
  }
  if (i == nRE)
    return;

  if (i < nRE-1) {
    a0[0] += x[0] * cf0;
    a1[0] += x[0] * cf1;
    x += incX;
    a0[1] += x[0] * cf0;
    a1[1] += x[0] * cf1;
    x += incX;
    a0 += 2;
    a1 += 2;
    i += 2;
  }
  if (i < nRE) {
    a0[0] += x[0] * cf0;
    a1[0] += x[0] * cf1;
    x += incX;
    i++;
  }
}

static void
_dmvec_vpur_rank(mdata_t *A, const mvec_t *X, const mvec_t *Y,  double alpha, 
                 int S, int L, int R, int E)
{
  register const double *y, *x;
  register double *Ac;
  int j, nRE, nSL;

  Ac = &A->md[S*A->step + R];
  y  = &Y->md[S*Y->inc];
  x  = &X->md[R*X->inc];
  nRE = E - R;
  nSL = L - S;

  for (j = 0; j < nSL-1; j += 2) {
    //Ar = Ac;
    _inner_mv_daxpy2(Ac, Ac+A->step, x, X->inc, y, y+Y->inc, alpha, nRE);
    y += 2*Y->inc;
    Ac += 2*A->step;
  }
  if (j == nSL)
    return;
  if (j < nSL) {
    _inner_mv_daxpy(Ac, x, X->inc, y, alpha, nRE);
  }
}

// A = A + alpha * x * y.T; A is M*N, x is M*1, Y is N*1
void dmvec_rank(mdata_t *A, const mvec_t *X,  const mvec_t *Y, double alpha, 
                int S, int L, int R, int E,
                int NB, int MB)
{
  int i, j, nI, nJ;

  if (NB <= 0) {
    NB = L - S;
  }
  if (MB <= 0) {
    MB = E - R;
  }
  for (j = S; j < L; j += NB) {
    nJ = L - j < NB ? L - j : NB;
    for (i = R; i < E; i += MB) {
      nI = E - i < MB ? E - i : MB;
      _dmvec_vpur_rank(A, X, Y, alpha, j, j+nJ, i, i+nI);
    }
  }
}


static void
_dmvec_vpur_syr_upper(mdata_t *A, const mvec_t *X0, const mvec_t *X1,
                      double alpha, int S, int L, int R, int E)
{
  // assert(S==R && L==E)
  register int i, j, nRE, nSL;
  register double *Ac, *Ar0, *Ar1, *At, cf;
  register const double *x0, *x1, *Xc;
  
  Xc = &X0->md[R*X0->inc];
  x1 = &X1->md[S*X1->inc];
  //At = &A->md[S*A->step+R];
  nRE = E - R;
  nSL = L - S;

  //printf("R=%d, E=%d, E-R: %d, S=%d, L=%d, L-S: %d\n", R, E, E-R, S, L, L-S);

  Ac = &A->md[S*A->step+R];
  for (j = 0; j < nSL-1; j += 2) {
    x0  = Xc;
    Ar0 = Ac;
    Ar1 = Ar0 + A->step;
    // add to 2 columns
    _inner_mv_daxpy2(Ac, Ac+A->step, x0, X0->inc, x1, x1+X1->inc, alpha, S+j+1);
    
    // diagonal entry on 2nd column is missing ... row j
    cf = alpha * x1[X1->inc];
    Ar1[S+j+1] += x0[(S+j+1)*X0->inc] * cf;
    
    //printf("A=\n"); print_tile(At, A->step, E-R, L-S);
    x1 += 2*X1->inc;
    Ac += 2*A->step;
  }
  if (j == nSL)
    return;
  if (j < nSL) {
    x0 = Xc;
    _inner_mv_daxpy(Ac, x0, X0->inc, x1, alpha, nRE);
  }
}

static void
_dmvec_vpur_syr_lower(mdata_t *A, const mvec_t *X0, const mvec_t *X1,
                      double alpha, int S, int L, int R, int E)
{
  // assert(S==R && L==E)
  register int i, j, nRE, nSL, nSE;
  register double *Ac, *Ar0, *Ar1, *At, cf;
  register const double *x0, *x1, *Xc;
  
  Xc = &X0->md[R*X0->inc];
  x1 = &X1->md[R*X1->inc];
  At = &A->md[S*A->step+R];
  nRE = E - R;
  nSL = L - S;
  nSE = E - S;

  //printf("start: Xc[%d] = %.1f\n", R, Xc[0]);
  Ac = &A->md[S*A->step + R];
  for (j = 0; j < nSL-1; j += 2) {
    x0 = Xc;
    Ar0 = Ac + j;
    //printf("j=%d, x0=%.1f, x1=%.1f\n", j, x0[0], x1[0]);

    // do diagonal entry on 1st column ... row j
    cf = alpha * x1[0];
    Ar0[0] += x0[0] * cf;
    Ar0++;
    x0 += X0->inc;

    // add to 2 columns
    _inner_mv_daxpy2(Ar0, Ar0+A->step, x0, X0->inc, x1, x1+X1->inc, alpha, nSE-j-1);
    
    //printf("A=\n"); print_tile(At, A->step, E-R, L-S);
    x1 += 2*X1->inc;
    Ac += 2*A->step;
    Xc += 2*X0->inc;
  }
  if (j == nSL)
    return;
  if (j < nSL) {
    x0 = Xc;
    //x1 = Xc;
    //printf("last column: j=%d, nSE=%d\n", j, nSE);
    _inner_mv_daxpy(Ac+j, x0, X0->inc, x1, alpha, nSE-j);
  }
}

// A = A + alpha * x * y.T; A is N*N symmetric lower|upper, x is N*1
void dmvec_symv_rank(mdata_t *A, const mvec_t *X,  double alpha, int flags,
                     int S, int L, int NB)
{
  int i, j, nI, nJ, sR, sE;

  if (L-S <= 0) {
    return;
  }
  if (NB <= 0) {
    NB = L - S;
  }
  if (flags & MTX_UPPER) {
    for (j = S; j < L; j += NB) {
      nJ = L - j < NB ? L - j : NB;
      _dmvec_vpur_syr_upper(A, X, X, alpha, j, j+nJ, 0, j+nJ);
    }
  } else {
    for (j = S; j < L; j += NB) {
      nJ = L - j < NB ? L - j : NB;
      _dmvec_vpur_syr_lower(A, X, X, alpha, j, j+nJ, j, L);
    }
  }
}

void dmvec_trmv_upd(mdata_t *A, const mvec_t *X,  const mvec_t *Y, double alpha, int flags,
                    int S, int L, int NB)
{
  int i, j, nI, nJ, sR, sE;

  if (L-S <= 0) {
    return;
  }
  if (NB <= 0) {
    NB = L - S;
  }
  if (flags & MTX_UPPER) {
    for (j = S; j < L; j += NB) {
      nJ = L - j < NB ? L - j : NB;
      _dmvec_vpur_syr_upper(A, X, Y, alpha, j, j+nJ, 0, j+nJ);
    }
  } else {
    for (j = S; j < L; j += NB) {
      nJ = L - j < NB ? L - j : NB;
      _dmvec_vpur_syr_lower(A, X, Y, alpha, j, j+nJ, j, L);
    }
  }
}

static void
_dmvec_vpur_syr2(mdata_t *A, const mvec_t *X, const mvec_t *Y,  double alpha, 
                 int flags, int S, int L, int R, int E)
{
  register int i, j, nRE, nSL;
  register double *Ac, *Ar, cf;
  register const double *y, *x;
  
  Ac = &A->md[S*A->step + R];
  y  = &Y->md[S*Y->inc];
  x  = &X->md[R*X->inc];
  nRE = E - R;
  nSL = L - S;

  for (j = 0; j < nSL-1; j += 2) {
    Ar = Ac;
    _inner_mv_daxpy2(Ac, Ac+A->step, x, X->inc, y, y+Y->inc, alpha, nRE);
    y += 2*Y->inc;
    Ac += 2*A->step;
  }
  if (j == nSL)
    return;
  if (j < nSL) {
    _inner_mv_daxpy(Ac, x, X->inc, y, alpha, nRE);
  }
}

// A = A + alpha * x * y.T; A is N*N symmetric lower|upper,
// x is N*1 or 1*N, Y is N*1 or 1*N
void dmvec_symv_rank2(mdata_t *A, const mvec_t *X,  const mvec_t *Y,
                      double alpha, int flags,
                      int S, int L, int NB)
{
  int i, j, nI, nJ, sR, sE;

  if (L - S <= 0) {
    return;
  }

  if (NB <= 0) {
    NB = L - S;
  }
  if (flags & MTX_UPPER) {
    for (j = S; j < L; j += NB) {
      nJ = L - j < NB ? L - j : NB;
      _dmvec_vpur_syr_upper(A, X, Y, alpha, j, j+nJ, 0, j+nJ);
      _dmvec_vpur_syr_upper(A, Y, X, alpha, j, j+nJ, 0, j+nJ);
    }
  } else {
    for (j = S; j < L; j += NB) {
      nJ = L - j < NB ? L - j : NB;
      _dmvec_vpur_syr_lower(A, X, Y, alpha, j, j+nJ, j, L);
      _dmvec_vpur_syr_lower(A, Y, X, alpha, j, j+nJ, j, L);
    }
  }
}

// Local Variables:
// indent-tabs-mode: nil
// End:

