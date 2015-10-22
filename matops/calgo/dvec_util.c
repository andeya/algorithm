
// Copyright (c) Harri Rautila, 2012,2013

// This file is part of github.com/hrautila/matops package. It is free software,
// distributed under the terms of GNU Lesser General Public License Version 3, or
// any later version. See the COPYING tile included in this archive.

#include <stdio.h>
#include <stdint.h>
#include <math.h>

#include "cmops.h"
#include "inner_vec_dot.h"
#include "inner_vec_axpy.h"

/*
 * Z[0] = alpha*(X*Y.T) + beta*Z[0]
 */
void dvec_dots(mvec_t *Z, const mvec_t *X,  const mvec_t *Y, double alpha, double beta, int N)
{
  register int i;
  register double c0, c1;
  register const double *Xc, *Yc;

  Xc = X->md;
  Yc = Y->md;
  c0 = 0.0; c1 = 0.0;
  for (i = 0; i < N-1; i += 2) {
    c0 += Xc[0] * Yc[0];
    Xc += X->inc;
    Yc += Y->inc;
    c1 += Xc[0] * Yc[0];
    Xc += X->inc;
    Yc += Y->inc;
  }    
  if (i < N) {
    c0 += Xc[0] * Yc[0];
  }
  Z->md[0] *= beta;
  Z->md[0] += alpha * (c0 + c1);
}

/*
 * return: alpha*(X*Y.T)
 */
double dvec_dot(const mvec_t *X,  const mvec_t *Y, double alpha, int N)
{
  register int i;
  register double c0, c1;
  register const double *Xc, *Yc;

  Xc = X->md;
  Yc = Y->md;
  c0 = 0.0; c1 = 0.0;
  for (i = 0; i < N-1; i += 2) {
    c0 += Xc[0] * Yc[0];
    Xc += X->inc;
    Yc += Y->inc;
    c1 += Xc[0] * Yc[0];
    Xc += X->inc;
    Yc += Y->inc;
  }    
  if (i < N) {
    c0 += Xc[0] * Yc[0];
  }
  return alpha * (c0 + c1);
}

/*
 * Y = alpha*X + Y
 */
void dvec_axpy(mvec_t *Y,  const mvec_t *X, double alpha, int N)
{
  register int i;
  register double *Yc;
  register const double *Xc;

  Xc = X->md;
  Yc = Y->md;
  for (i = 0; i < N-1; i += 2) {
    Yc[0]      += alpha * Xc[0];
    Yc[Y->inc] += alpha * Xc[X->inc];
    Xc += (X->inc << 1);
    Yc += (Y->inc << 1);
  }    
  if (i < N) {
    Yc[0] += alpha * Xc[0];
  }
}

/*
 * return ||X - Y||_2
 */
double dvec_diff_nrm2(const mvec_t *X,  const mvec_t *Y, int N)
{
  register int i;
  register double c0, c1, d0, d1;
  register const double *Xc, *Yc;

  Xc = X->md;
  Yc = Y->md;
  c0 = 0.0; c1 = 0.0;
  for (i = 0; i < N-1; i += 2) {
    d0 = fabs(Xc[0] - Yc[0]);
    c0 += d0 * d0;

    d1 = fabs(Xc[X->inc] - Yc[Y->inc]);
    c1 += d1 * d1;
    Xc += (X->inc << 1);
    Yc += (Y->inc << 1);
  }    
  if (i < N) {
    d0 = fabs(Xc[0] - Yc[0]);
    c0 += d0 * d0;
  }
  return sqrt(c0 + c1);
}

// return vector norm 
double dvec_nrm2(const mvec_t *X,  int N)
{
  register int i;
  register double c0, c1, d0, d1, absx;
  register const double *Xc;

  Xc = X->md;
  c0 = 0.0; c1 = 0.0;
  for (i = 0; i < N-1; i += 2) {
    absx = fabs(Xc[0]);
    c0 += absx * absx;

    absx = fabs(Xc[X->inc]);
    c1 += absx * absx;
    Xc += (X->inc << 1);
  }    
  if (i < N) {
    absx = fabs(Xc[0]);
    c0 += absx * absx;
  }
  return sqrt(c0 + c1);
}

// return sum of absolute values
double dvec_asum(const mvec_t *X,  int N)
{
  register int i;
  register double c0, c1;
  register const double *Xc;

  Xc = X->md;
  c0 = 0.0; c1 = 0.0;
  for (i = 0; i < N-1; i += 2) {
    c0 += fabs(Xc[0]);
    c1 += fabs(Xc[X->inc]);
    Xc += (X->inc << 1);
  }    
  if (i < N) {
    c0 += fabs(Xc[0]);
  }
  return c0 + c1;
}

// return index of max absolute value
int dvec_iamax(const mvec_t *X,  int N)
{
  register int i, ix, n;
  register double max, c0, c1;
  register const double *Xc;

  if (N <= 1)
    return 0;
  Xc = X->md;
  max = 0.0;
  ix = 0;
  //N--;
  for (i = 0; i < N-1; i += 2) {
    c0 = fabs(Xc[0]);
    c1 = fabs(Xc[X->inc]);
    if (c1 > c0) {
      n = 1;
      c0 = c1;
    }
    if (c0 > max) {
      ix = i+n;
      max = c0;
    }
    Xc += (X->inc << 1);
    n = 0;
  }    
  if (i < N) {
    c0 = fabs(Xc[0]);
    ix = c0 > max ? N-1 : ix;
  }
  return ix;
}

/*
 * X <--> Y
 */
void dvec_swap(mvec_t *X,  mvec_t *Y, int N)
{
  register int i;
  register double tmp;
  register double *Xc, *Yc;

  Xc = X->md;
  Yc = Y->md;
  for (i = 0; i < N-1; i += 2) {
    tmp = Xc[0];
    Xc[0] = Yc[0];
    Yc[0] = tmp;
    Xc += X->inc;
    Yc += Y->inc;
    tmp = Xc[0];
    Xc[0] = Yc[0];
    Yc[0] = tmp;
    Xc += X->inc;
    Yc += Y->inc;
  }    
  if (i < N) {
    tmp = Xc[0];
    Xc[0] = Yc[0];
    Yc[0] = tmp;
  }
}

/*
 * X = X/alpha
 */
void dvec_invscal(mvec_t *X,  double alpha, int N)
{
  register int i;
  register double *Xc;

  Xc = X->md;
  for (i = 0; i < N-1; i += 2) {
    Xc[0] /= alpha;
    Xc += X->inc;
    Xc[0] /= alpha;
    Xc += X->inc;
  }    
  if (i < N) {
    Xc[0] /= alpha;
  }
}

/*
 * X = alpha*X
 */
void dvec_scal(mvec_t *X,  double alpha, int N)
{
  register int i;
  register double *Xc;

  Xc = X->md;
  for (i = 0; i < N-1; i += 2) {
    Xc[0] *= alpha;
    Xc += X->inc;
    Xc[0] *= alpha;
    Xc += X->inc;
  }    
  if (i != N) {
    Xc[0] *= alpha;
  }
}

/*
 * X := Y
 */
void dvec_copy(mvec_t *X,  mvec_t *Y, int N)
{
  register int i;
  register double *Xc, *Yc;

  Xc = X->md;
  Yc = Y->md;
  for (i = 0; i < N-1; i += 2) {
    Xc[0] = Yc[0];
    Xc[X->inc] = Yc[Y->inc];
    Xc += 2*X->inc;
    Yc += 2*Y->inc;
  }    
  if (i < N) {
    Xc[0] = Yc[0];
  }
}

// Local Variables:
// indent-tabs-mode: nil
// End:
