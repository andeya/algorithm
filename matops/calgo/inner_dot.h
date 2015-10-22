
// Copyright (c) Harri Rautila, 2012,2013

// This file is part of github.com/hrautila/matops package. It is free software,
// distributed under the terms of GNU Lesser General Public License Version 3, or
// any later version. See the COPYING tile included in this archive.

#ifndef __INNER_DOT_SSENEW
#define __INNER_DOT_SSENEW

#include <x86intrin.h>


static inline 
void _inner_ddot_ssen(double *Cr, const double *Ar, const double *Br, double alpha, int nVP)
{
  register int k;
  register double cval = 0.0;
  register __m128d A0, B0, C0, F0, ALP;

  C0 = _mm_set1_pd(0.0);
  ALP = _mm_set1_pd(alpha);

  // unrolling of loops;
  for (k = nVP-3; k > 0; k -= 4) {
    A0 = _mm_load_pd(Ar);
    B0 = _mm_load_pd(Br);
    F0 = A0 * B0;
    C0 = C0 + F0;
    A0 = _mm_load_pd(&Ar[2]);
    B0 = _mm_load_pd(&Br[2]);
    F0 = A0 * B0;
    C0 = C0 + F0;
    Ar += 4;
    Br += 4;
  }
  k += 3;
  switch (k) {
  case 3:
  case 2:
    A0 = _mm_load_pd(Ar);
    B0 = _mm_load_pd(Br);
    F0 = A0 * B0;
    C0 = C0 + F0;
    Ar += 2; Br += 2;
    if (k == 2)
      goto update;
  case 1:
    cval = Ar[0] * Br[0];
    Cr[0] += cval * alpha;
  }
 update:
  C0 = C0 * ALP;
  Cr[0] += C0[1];
  Cr[0] += C0[0];
}


static inline 
void _inner_ddot2_ssen(double *Cr0, double *Cr1, const double *Ar,
		       const double *b0, const double *b1, double alpha, int nVP)
{
  register int k, i;
  register double cval0, cval1;
  register __m128d A0, B0, C0, F0, A1, B1, C1, F1, ALP;

  C0 = _mm_set1_pd(0.0);
  C1 = _mm_set1_pd(0.0);
  ALP = _mm_set1_pd(alpha);

  // unrolling of loops;
  for (k = nVP-1; k > 0; k -= 2) {
    A0 = _mm_load_pd(Ar);
    B0 = _mm_load_pd(b0);
    B1 = _mm_load_pd(b1);
    F0 = A0 * B0;
    C0 = C0 + F0;
    F1 = A0 * B1;
    C1 = C1 + F1;
    Ar += 2;
    b0 += 2;
    b1 += 2;
  }
  if (k == 0) {
    cval0 = Ar[0] * b0[0];
    cval1 = Ar[0] * b1[0];
    Cr0[0] += cval0 * alpha;
    Cr1[0] += cval1 * alpha;
  }
  // SSE3 horizontal add
  F0 = _mm_hadd_pd(C0, C1);
  Cr0[0] += F0[1] * alpha;
  Cr1[0] += F0[0] * alpha;
  //C0 = C0 * ALP;
  //C1 = C1 * ALP;
  //Cr0[0] += C0[0];
  //Cr0[0] += C0[1];
  //Cr1[0] += C1[0];
  //Cr1[0] += C1[1];
}


// update Cr0 += sum(Ar * b0); Cr1 += sum(Ar * b1); ...
static inline 
void _inner_ddot4_ssen(double *Cr0, double *Cr1, double *Cr2, double *Cr3,
		       const double *Ar, const double *b0, const double *b1,
		       const double *b2, const double *b3, double alpha, int nVP)
{
  register int k, i;
  register double cval0, cval1, cval2, cval3;
  register __m128d A0, B0, C0, F0, A1, B1, C1, F1;
  register __m128d A2, B2, C2, F2, A3, B3, C3, F3, ALP;

  C0 = _mm_set1_pd(0.0);
  C1 = C0;
  C2 = C0;
  C3 = C0;
  ALP = _mm_set1_pd(alpha);

  // unrolling of loops;
  for (k = nVP-1; k > 0; k -= 2) {
    A0 = _mm_load_pd(Ar);
    B0 = _mm_load_pd(b0);
    B1 = _mm_load_pd(b1);
    B2 = _mm_load_pd(b2);
    B3 = _mm_load_pd(b3);
    F0 = A0 * B0;
    C0 = C0 + F0;
    F1 = A0 * B1;
    C1 = C1 + F1;
    F2 = A0 * B2;
    C2 = C2 + F2;
    F3 = A0 * B3;
    C3 = C3 + F3;

    Ar += 2;
    b0 += 2;
    b1 += 2;
    b2 += 2;
    b3 += 2;
  }
  if (k == 0) {
    cval0 = Ar[0] * b0[0];
    cval1 = Ar[0] * b1[0]; 
    cval2 = Ar[0] * b2[0];
    cval3 = Ar[0] * b3[0]; 
    Cr0[0] += cval0 * alpha;
    Cr1[0] += cval1 * alpha;
    Cr2[0] += cval2 * alpha;
    Cr3[0] += cval3 * alpha;
  }
 update:
  F0 = _mm_hadd_pd(C0, C1);
  F1 = _mm_hadd_pd(C2, C3);
  Cr0[0] += F0[1] * alpha;
  Cr1[0] += F0[0] * alpha;
  Cr2[0] += F1[1] * alpha;
  Cr3[0] += F1[0] * alpha;

  //C0 = C0 * ALP;
  //C1 = C1 * ALP;
  //C2 = C2 * ALP;
  //C3 = C3 * ALP;
  //Cr0[0] += C0[0];
  //Cr0[0] += C0[1];
  //Cr1[0] += C1[0];
  //  Cr1[0] += C1[1];
  //Cr2[0] += C2[0];
  //Cr2[0] += C2[1];
  //Cr3[0] += C3[0];
  //Cr3[0] += C3[1];
}


// update Cr0[0] += sum(a0 * b0); Cr0[1] += sum(a1 * b0); Cr1[0] += sum(a0 * b1); ....
static inline 
void _inner_ddot4_2_ssen(double *Cr0, double *Cr1, double *Cr2, double *Cr3,
			 const double *a0, const double *a1,
			 const double *b0, const double *b1,
			 const double *b2, const double *b3, double alpha, int nVP)
{
  register int k, i;
  register double cval0, cval1, cval2, cval3;
  register __m128d A0, B0, C0_0, F0, A1, B1, C1_0, F1;
  register __m128d A2, B2, C2_0, F2, A3, B3, C3_0, F3, ALP;
  register __m128d C0_1, C1_1, C2_1, C3_1;

  C0_0 = _mm_set1_pd(0.0);
  C1_0 = _mm_set1_pd(0.0);
  C2_0 = _mm_set1_pd(0.0);
  C3_0 = _mm_set1_pd(0.0);
  C0_1 = C0_0;
  C1_1 = C1_0;
  C2_1 = C2_0;
  C3_1 = C3_0;
  ALP = _mm_set1_pd(alpha);

  // unrolling of loops;
  for (k = nVP-1; k > 0; k -= 2) {
    A0 = _mm_load_pd(a0);
    B0 = _mm_load_pd(b0);
    B1 = _mm_load_pd(b1);
    B2 = _mm_load_pd(b2);
    B3 = _mm_load_pd(b3);
    F0   = A0 * B0;
    C0_0 = C0_0 + F0;
    F1   = A0 * B1;
    C1_0 = C1_0 + F1;
    F2   = A0 * B2;
    C2_0 = C2_0 + F2;
    F3   = A0 * B3;
    C3_0 = C3_0 + F3;

    A0 = _mm_load_pd(a1);
    F0   = A0 * B0;
    C0_1 = C0_1 + F0;
    F1   = A0 * B1;
    C1_1 = C1_1 + F1;
    F2   = A0 * B2;
    C2_1 = C2_1 + F2;
    F3   = A0 * B3;
    C3_1 = C3_1 + F3;

    a0 += 2;
    a1 += 2;
    b0 += 2;
    b1 += 2;
    b2 += 2;
    b3 += 2;
  }
  if (k == 0) {
    cval0 = a0[0] * b0[0];
    cval1 = a0[0] * b1[0]; 
    cval2 = a0[0] * b2[0];
    cval3 = a0[0] * b3[0]; 
    Cr0[0] += cval0 * alpha;
    Cr1[0] += cval1 * alpha;
    Cr2[0] += cval2 * alpha;
    Cr3[0] += cval3 * alpha;
    cval0 = a1[0] * b0[0];
    cval1 = a1[0] * b1[0]; 
    cval2 = a1[0] * b2[0];
    cval3 = a1[0] * b3[0]; 
    Cr0[1] += cval0 * alpha;
    Cr1[1] += cval1 * alpha;
    Cr2[1] += cval2 * alpha;
    Cr3[1] += cval3 * alpha;
  }
 update:
  F0 = _mm_hadd_pd(C0_0, C1_0);
  F1 = _mm_hadd_pd(C2_0, C3_0);
  Cr0[0] += F0[1] * alpha;
  Cr1[0] += F0[0] * alpha;
  Cr2[0] += F1[1] * alpha;
  Cr3[0] += F1[0] * alpha;

  F2 = _mm_hadd_pd(C0_1, C1_1);
  F3 = _mm_hadd_pd(C2_1, C3_1);
  Cr0[1] += F2[1] * alpha;
  Cr1[1] += F2[0] * alpha;
  Cr2[1] += F3[1] * alpha;
  Cr3[1] += F3[0] * alpha;

  /*
  C0_0 = C0_0 * ALP;
  C1_0 = C1_0 * ALP;
  C2_0 = C2_0 * ALP;
  C3_0 = C3_0 * ALP;
  Cr0[0] += C0_0[0];
  Cr0[0] += C0_0[1];
  Cr1[0] += C1_0[0];
  Cr1[0] += C1_0[1];
  Cr2[0] += C2_0[0];
  Cr2[0] += C2_0[1];
  Cr3[0] += C3_0[0];
  Cr3[0] += C3_0[1];

  C0_1 = C0_1 * ALP;
  C1_1 = C1_1 * ALP;
  C2_1 = C2_1 * ALP;
  C3_1 = C3_1 * ALP;
  Cr0[1] += C0_1[0];
  Cr0[1] += C0_1[1];
  Cr1[1] += C1_1[0];
  Cr1[1] += C1_1[1];
  Cr2[1] += C2_1[0];
  Cr2[1] += C2_1[1];
  Cr3[1] += C3_1[0];
  Cr3[1] += C3_1[1];
  */
}


#endif

// Local Variables:
// indent-tabs-mode: nil
// End:
