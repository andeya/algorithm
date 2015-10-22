
// Copyright (c) Harri Rautila, 2012,2013

// This file is part of github.com/hrautila/matops package. It is free software,
// distributed under the terms of GNU Lesser General Public License Version 3, or
// any later version. See the COPYING tile included in this archive.

#ifndef __INNER_VEC_DOT_H
#define __INNER_VEC_DOT_H 1

#include <x86intrin.h>
#include <emmintrin.h>

// Update Y one cell with A[i:]*X
static inline
void _inner_vec_ddot(double *y0, int incY, const double *a0,
                     const double *x0, int incX, double alpha, int nC)
{
  register int i;
  register double ytmp;

  ytmp = 0.0;
  for (i = 0; i < nC-3; i += 4) {
    ytmp += a0[0] * x0[0];
    //x0 += incX;
    ytmp += a0[1] * x0[incX];
    //x0 += incX;
    ytmp += a0[2] * x0[2*incX];
    //x0 += incX;
    ytmp += a0[3] * x0[3*incX];
    //x0 += incX;
    x0 += incX << 2;
    a0 += 4;
  }
  if (i == nC) 
    goto update;

  if (i < nC-1) {
    ytmp += a0[0] * x0[0];
    //x0 += incX;
    ytmp += a0[1] * x0[incX];
    //x0 += incX;
    x0 += incX << 1;
    a0 += 2;
    i += 2;
  }
  if (i < nC) {
    ytmp += a0[0] * x0[0];
    x0 += incX;
    i++;
  }
 update:
  y0[0] += ytmp * alpha;
  
}

static inline
void _inner_vec2_ddot(double *y0, int incY, const double *a0, const double *a1,
                      const double *x0, int incX, double alpha, int nC)
{
  register int i;
  register double ytmp0, ytmp1;

  ytmp0 = 0.0;
  ytmp1 = 0.0;
  for (i = 0; i < nC-3; i += 4) {
    ytmp0 += a0[0] * x0[0];
    ytmp1 += a1[0] * x0[0];
    //x0 += incX;
    ytmp0 += a0[1] * x0[incX];
    ytmp1 += a1[1] * x0[incX];
    //x0 += incX;
    ytmp0 += a0[2] * x0[2*incX];
    ytmp1 += a1[2] * x0[2*incX];
    //x0 += incX;
    ytmp0 += a0[3] * x0[3*incX];
    ytmp1 += a1[3] * x0[3*incX];
    //x0 += incX;
    x0 += incX << 2;
    a0 += 4;
    a1 += 4;
  }
  if (i == nC) 
    goto update;

  if (i < nC-1) {
    ytmp0 += a0[0] * x0[0];
    ytmp1 += a1[0] * x0[0];
    //x0 += incX;
    ytmp0 += a0[1] * x0[incX];
    ytmp1 += a1[1] * x0[incX];
    //x0 += incX;
    x0 += incX << 1;
    a0 += 2;
    a1 += 2;
    i += 2;
  }
  if (i < nC) {
    ytmp0 += a0[0] * x0[0];
    ytmp1 += a1[0] * x0[0];
    x0 += incX;
    i++;
  }
 update:
  y0[0] += ytmp0 * alpha;
  y0[incY] += ytmp1 * alpha;
  
}

static inline
void _inner_vec_ddot_sse(double *y0, int incY, const double *a0,
                         const double *x0, double alpha, int nC, int oddstart)
{
  register int i;
  register double ytmp;
  register __m128d Y0, Y1, A0, X0, TMP0, TMP1;

  ytmp = 0.0;
  Y0 = _mm_set1_pd(0.0);
  Y1 = Y0;

  if (oddstart) {
    ytmp += a0[0] * x0[0];
    x0++;
    a0++;
    nC--;
  }

  for (i = 0; i < nC-3; i += 4) {
    A0 = _mm_load_pd(a0);
    X0 = _mm_load_pd(x0);
    TMP0 = A0 * X0;
    Y0 = Y0 + TMP0;
    x0 += 2;
    a0 += 2;

    A0 = _mm_load_pd(a0);
    X0 = _mm_load_pd(x0);
    TMP1 = A0 * X0;
    Y1 = Y1 + TMP1;
    x0 += 2;
    a0 += 2;
  }
  if (i == nC) 
    goto update;

  if (i < nC-1) {
    A0 = _mm_load_pd(a0);
    X0 = _mm_load_pd(x0);
    TMP0 = A0 * X0;
    Y0 = Y0 + TMP0;
    x0 += 2;
    a0 += 2;
    i += 2;
  }
  if (i < nC) {
    ytmp += a0[0] * x0[0];
    i++;
  }
 update:
  Y0 = _mm_hadd_pd(Y0, Y1);
  ytmp += Y0[0];
  ytmp += Y0[1];
  ytmp *= alpha;
  y0[0] += ytmp;
  
}

static inline
void _inner_vec2_ddot_sse(double *y0, int incY, const double *a0, const double *a1,
                         const double *x0, double alpha, int nC, int oddstart)
{
  register int i;
  register double ytmp0, ytmp1;
  register __m128d Y0, Y1, A0, A1, X0, TMP0, TMP1;

  ytmp0 = 0.0;
  ytmp1 = 0.0;
  Y0 = _mm_set1_pd(0.0);
  Y1 = _mm_set1_pd(0.0);

  if (oddstart) {
    ytmp0 += a0[0] * x0[0];
    ytmp1 += a1[0] * x0[0];
    x0++;
    a0++;
    a1++;
    nC--;
  }

  for (i = 0; i < nC-3; i += 4) {
    A0 = _mm_load_pd(a0);
    A1 = _mm_load_pd(a1);
    X0 = _mm_load_pd(x0);
    TMP0 = A0 * X0;
    TMP1 = A1 * X0;
    Y0 = Y0 + TMP0;
    Y1 = Y1 + TMP1;
    x0 += 2;
    a0 += 2;
    a1 += 2;

    A0 = _mm_load_pd(a0);
    A1 = _mm_load_pd(a1);
    X0 = _mm_load_pd(x0);
    TMP0 = A0 * X0;
    TMP1 = A1 * X0;
    Y0 = Y0 + TMP0;
    Y1 = Y1 + TMP1;
    x0 += 2;
    a0 += 2;
    a1 += 2;
  }
  if (i == nC) 
    goto update;

  if (i < nC-1) {
    A0 = _mm_load_pd(a0);
    A1 = _mm_load_pd(a1);
    X0 = _mm_load_pd(x0);
    TMP0 = A0 * X0;
    TMP1 = A1 * X0;
    Y0 = Y0 + TMP0;
    Y1 = Y1 + TMP1;
    x0 += 2;
    a0 += 2;
    a1 += 2;
    i += 2;
  }
  if (i < nC) {
    ytmp0 += a0[0] * x0[0];
    ytmp1 += a1[0] * x0[0];
    i++;
  }
 update:
  TMP1 = _mm_hadd_pd(Y0, Y1);
  ytmp0 += TMP1[0];
  ytmp1 += TMP1[1];
  y0[0]    += ytmp0 * alpha;
  y0[incY] += ytmp1 * alpha;
  
}

#endif


// Local Variables:
// indent-tabs-mode: nil
// End:
