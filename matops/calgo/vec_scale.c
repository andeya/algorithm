
// Copyright (c) Harri Rautila, 2012,2013

// This file is part of github.com/hrautila/matops package. It is free software,
// distributed under the terms of GNU Lesser General Public License Version 3, or
// any later version. See the COPYING tile included in this archive.

#include <stdio.h>
#include "cmops.h"

// Kahan summation for DOT product:
//    http://en.wikipedia.org/wiki/Kahan_summation_algorithm
// Not usually used in BLAS because of perfomance considerations.

// this is __not__ kahan summation 
double ddot_vec(const double *X, const double *Y, int incX, int incY, int N)
{
  register int i;
  register double c0, c1;
  c0 = 0.0; c1 = 0.0;
  for (i = 0; i < N-1; i += 2) {
    c0 += X[0] * Y[0];
    X += incX;
    Y += incY;
    c1 += X[0] * Y[0];
    X += incX;
    Y += incY;
  }    
  if (i == N) {
    return c0 + c1;
  }
  c0 += X[0] * Y[0];
  return c0 + c1;
}

// Y = alpha * X + beta * Y;
void daxpy_vec(double *Y, const double *X, double alpha, double beta,
               int incX, int incY, int N)
{
  // could make SSE version if incX == incY == 1 and align(Y) == align(X)

  register int i;
  register double *y0;
  register const double *x0;

  if (beta != 1.0) {
    y0 = Y;
    if (beta == 0.0) {
      for (i = 0; i < N; i++) {
        y0[0] = 0.0;
        y0++;
      }
    } else {
      for (i = 0; i < N; i++) {
        y0[0] *= beta;
        y0++;
      }
    }
  }
  y0 = Y; x0 = X;
  for (i = 0; i < N-3; i += 4) {
    y0[0] += alpha * x0[0];
    x0 += incX;
    y0 += incY;
    y0[0] += alpha * x0[0];
    x0 += incX;
    y0 += incY;
    y0[0] += alpha * x0[0];
    x0 += incX;
    y0 += incY;
    y0[0] += alpha * x0[0];
    x0 += incX;
    y0 += incY;
  }    
  if (i == N) {
    return;
  }
  if (i < N-1) {
    y0[0] += alpha * x0[0];
    x0 += incX;
    y0 += incY;
    y0[0] += alpha * x0[0];
    x0 += incX;
    y0 += incY;
    i += 2;
  }
  if (i < N) {
    y0[0] += alpha * x0[0];
  }
  return;
}

// Scale a vector of N elements with incX interval.
void dscale_vec(double *X, int incX, double f0, int N)
{
  register int i;
  if (f0 == 1.0) {
    return;
  }
  if (f0 == 0.0) {
    for (i = 0; i < N; i++) {
      X[0] = 0.0;
      X += incX;
    }
    return;
  }
  for (i = 0; i < N; i++) {
    X[0] *= f0;
    X += incX;
  }
}

// Scale a tile of M rows by N columns with leading index ldX.
void dscale_tile(double *X, int ldX, double f0, int M, int N)
{
  register double *Xr, *Xc;
  register int i, j;
  if (f0 == 1.0) {
    return;
  }
  Xc = X;
  // set to zero
  if (f0 == 0.0) {
    for (j = 0; j < N; j++) {
      Xr = Xc;
      for (i = 0; i < M; i++) {
        Xr[0] = 0.0;
        Xr++;
      }
      Xc += ldX;
    }
    return;
  }

  // scale here
  for (j = 0; j < N; j++) {
    Xr = Xc;
    for (i = 0; i < M-3; i += 4) {
      Xr[0] *= f0;
      Xr[1] *= f0;
      Xr[2] *= f0;
      Xr[3] *= f0;
      Xr += 4;
    }
    if (i == M)
      goto increment;
    if (i < M-1) {
      Xr[0] *= f0;
      Xr[1] *= f0;
      Xr += 2;
      i += 2;
    }
    if (i < M) {
      Xr[0] *= f0;
      i++;
    }
  increment:
    Xc += ldX;
  }
  return;

}

// scale triangular; if MTX_UNIT set does not scale diagonal element
void dscale_triul(double *X, int ldX, double f0, int N, int flags)
{
  register double *Xr, *Xc;
  register int i, j;
  register int unit = flags & MTX_UNIT ? 1 : 0;
  if (f0 == 1.0) {
    return;
  }
  Xc = X;
  // set to zero
  if (f0 == 0.0) {
    if (flags & MTX_LOWER) {
      for (j = 0; j < N; j++) {
        Xr = Xc + j + unit;
        for (i = N; i > j+unit; i--) {
          Xr[0] = 0.0;
          Xr++;
        }
        Xc += ldX;
      }
    } else {
      for (j = 0; j < N; j++) {
        Xr = Xc;
        for (i = 0; i < j+unit; i++) {
          Xr[0] = 0.0;
          Xr++;
        }
        Xc += ldX;
      }
    }
    return;
  }

  // scale here
  if (flags & MTX_LOWER) {
    for (j = 0; j < N; j++) {
      Xr = Xc + j + unit;
      for (i = N; i > j+unit; i--) {
        Xr[0] *= f0;
        Xr ++;
      }
      Xc += ldX;
    }
  } else {
    for (j = 0; j < N; j++) {
      Xr = Xc;
      for (i = 0; i < j+unit; i ++) {
        Xr[0] *= f0;
        Xr ++;
      }
      Xc += ldX;
    }
  }
  return;

}

void print_tile(const double *D, int ldD, int nR, int nC)
{
  register int i, j;
  for (i = 0; i < nR; i++) {
    printf("[");
    for (j = 0; j < nC; j++) {
      if (j != 0)
        printf(", ");
      printf("%9.2e", D[j*ldD+i]);
    }
    printf("]\n");
  }
}

// Local Variables:
// indent-tabs-mode: nil
// End:

