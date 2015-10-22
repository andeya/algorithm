
// Copyright (c) Harri Rautila, 2012,2013

// This file is part of github.com/hrautila/matops package. It is free software,
// distributed under the terms of GNU Lesser General Public License Version 3, or
// any later version. See the COPYING tile included in this archive.


#ifndef _COLCPY_H
#define _COLCPY_H

extern void *memcpy(void *, const void *, size_t);

static inline
void colcpy(double * __restrict dst, int ldD,
            const double * __restrict src, int ldS, int nR, int nC)
{
  register int i;
  for (i = 0; i < nC; i++) {
    memcpy(dst, src, nR*sizeof(double));
    dst += ldD;
    src += ldS;
  }
}

// Copy nC columns of length nR from source to dest.
// Source row stride is ldS and destination row stride is ldD.
static inline
void colcpy4(double * __restrict dst, int ldD,
             const double * __restrict src, int ldS, int nR, int nC)
{
  register int i;
  for (i = 0; i < nC-3; i += 4) {
    memcpy(dst, src, nR*sizeof(double));
    dst += ldD;
    src += ldS;
    memcpy(dst, src, nR*sizeof(double));
    dst += ldD;
    src += ldS;
    memcpy(dst, src, nR*sizeof(double));
    dst += ldD;
    src += ldS;
    memcpy(dst, src, nR*sizeof(double));
    dst += ldD;
    src += ldS;
  }
  if (i == nC)
    return;
  if (i < nC-1) {
    memcpy(dst, src, nR*sizeof(double));
    dst += ldD;
    src += ldS;
    memcpy(dst, src, nR*sizeof(double));
    dst += ldD;
    src += ldS;
    i += 2;
  }
  if (i < nC) {
    memcpy(dst, src, nR*sizeof(double));
  }
}


static inline
void colcpy_trans(double * __restrict dst, int ldD,
                  const double * __restrict src, int ldS, int nR, int nC)
{
  register double *Dc, *Dr;
  register const double *Sc, *Sr;
  register int j, i;
  Dc = dst; Sc = src;
  for (j = 0; j < nC; j++) {
    Dr = Dc;
    Sr = Sc;
    // incrementing Dr with ldD follows the dst row
    // and incrementing Sr with one follows the column
    for (i = 0; i <nR; i++) {
      *Dr = *Sr;
      Dr += ldD;
      Sr++;
    }
    // moves Dc pointer to next row on dst
    Dc++;
    // moves Sc pointer to next column on src
    Sc += ldS;
  }
}

// Copying 2 columns parallel distributes cache misses more evenly.
static inline
void colcpy2_trans(double * __restrict dst, int ldD,
		  const double * __restrict src, int ldS, int nR, int nC)
{
  register double *Dc, *Dr0, *Dr1;
  register const double *Sc, *Sr0, *Sr1;
  register int j, i;
  Dc = dst; Sc = src;
  for (j = 0; j < nC-1; j += 2) {
    Dr0 = Dc;
    Dr1 = Dr0 + 1;
    Sr0 = Sc;
    Sr1 = Sr0 + ldS;
    // incrementing Dr with ldD follows the dst row
    // and incrementing Sr with one follows the column
    for (i = nR; i > 0; i--) {
      Dr0[0] = Sr0[0];
      Dr1[0] = Sr1[0];
      Dr0 += ldD;
      Dr1 += ldD;
      Sr0 += 1;
      Sr1 += 1;
    }
    // moves Dc pointer to next row on dst
    Dc += 2;
    // moves Sc pointer to next column on src
    Sc += (ldS << 1);
  }
  if (j < nC) {
    Sr0 = Sc;
    Dr0 = Dc;
    for (i = nR; i > 0; i--) {
      Dr0[0] = Sr0[0];
      Dr0 += ldD;
      Sr0 += 1;
    }
  }
}

// Copying 4 columns parallel distributes cache misses even more evenly.
static inline
void colcpy4_trans(double * __restrict dst, int ldD,
		  const double * __restrict src, int ldS, int nR, int nC)
{
  register double *Dc, *dr0, *dr1, *dr2, *dr3;
  register const double *Sc, *sr0, *sr1, *sr2, *sr3;
  register int j, i;
  Dc = dst; Sc = src;
  for (j = 0; j < nC-3; j += 4) {
    dr0 = Dc;
    dr1 = dr0 + 1;
    dr2 = dr1 + 1;
    dr3 = dr2 + 1;
    sr0 = Sc;
    sr1 = sr0 + ldS;
    sr2 = sr1 + ldS;
    sr3 = sr2 + ldS;
    // incrementing Dr with ldD follows the dst row
    // and incrementing Sr with one follows the column
    for (i = nR; i > 0; i--) {
      dr0[0] = sr0[0];
      dr1[0] = sr1[0];
      dr2[0] = sr2[0];
      dr3[0] = sr3[0];
      dr0 += ldD;
      dr1 += ldD;
      dr2 += ldD;
      dr3 += ldD;
      sr0 += 1;
      sr1 += 1;
      sr2 += 1;
      sr3 += 1;
    }
    // moves Dc pointer to next row on dst
    Dc += 4;
    // moves Sc pointer to next column on src
    Sc += (ldS << 2);
  }
  if (j < nC-1) {
    dr0 = Dc;
    dr1 = dr0 + 1;
    sr0 = Sc;
    sr1 = sr0 + ldS;
    for (i = nR; i > 0; i--) {
      dr0[0] = sr0[0];
      dr1[0] = sr1[0];
      dr0 += ldD;
      dr1 += ldD;
      sr0 += 1;
      sr1 += 1;
    }
    Dc += 2;
    Sc += (ldS << 1);
    j += 2;
  }
  if (j < nC) {
    sr0 = Sc;
    dr0 = Dc;
    for (i = nR; i > 0; i--) {
      dr0[0] = sr0[0];
      dr0 += ldD;
      sr0 += 1;
    }
  }
}


// Copy upper tridiagonal and fill lower part to form full symmetric matrix
// result is symmetric matrix A and A = A.T
static inline
void colcpy_fill_low(double *dst, int ldD, const double *src, int ldS, int nR, int nC, int unit)
{
  //assert(nR == nC);
  register double *Dcu, *Dcl, *Drl, *Dru;
  register const double *Sc, *Sr;
  register int j, i;
  Dcu = dst; Sc = src;
  Dcl = dst;
  // fill dst row and column at the same time, following src columns
  for (j = 0; j < nC; j++) {
    Dru = Dcu;
    Drl = Dcl;
    Sr = Sc;
    for (i = 0; i < j; i++) {
      *Dru = *Sr;
      *Drl = *Sr;
      Sr++;
      Dru++; 
      Drl += ldD;
    }
    *Dru = unit ? 1.0 : *Sr;  // copy the diagonal entry
    Sc += ldS;                // next column in source
    Dcu += ldD;               // next column for upper triagonal
    Dcl++;                    // next row for lower triagonal
  }
}

// Copy lower tridiagonal and fill upper part to form full symmetric matrix;
// result is symmetric matrix A and A = A.T
static inline
void colcpy_fill_up(double *dst, int ldD, const double *src, int ldS, int nR, int nC, int unit)
{
  //assert(nR == nC);
  register double *Dcu, *Dcl, *Drl, *Dru;
  register const double *Sc, *Sr;
  register int j, i;
  Dcu = dst; Sc = src;
  Dcl = dst;
  // fill dst row and column at the same time, following src columns
  for (j = 0; j < nC; j++) {
    // start at same point and diverge down (Drl) and right (Dru)
    Dru = Dcu + j;
    Drl = Dcl + j;
    // start of data on column, j'th row (diagonal entry)
    Sr = Sc + j;
    // diagonal entry
    *Dru = unit ? 1.0 : *Sr;
    Sr++;
    Dru += ldD;
    Drl++;     
    // off diagonal entries
    for (i = 1; i < nC-j; i++) {
      *Dru = *Sr;
      *Drl = *Sr;
      Sr++;
      Dru += ldD;       // next column in row
      Drl++;            // next row in column 
    }
    // NEXT column in source
    Sc += ldS;
    // next column for upper triagonal
    Dcu += ldD;
    // next column for lower triagonal
    Dcl += ldD;
  }
}




// Transpose upper tridiagonal to fill lower part and fill dst upper part with zeros.
// Result is lower tridiagonal matrix with zero upper part.
static inline
void colcpy_trans_upper_and_zero(double *dst, int ldD, const double *src,
                                 int ldS, int nR, int nC, int unit)
{
  //assert(nR == nC);
  register double *Dcu, *Dcl, *Drl, *Dru;
  register const double *Sc, *Sr;
  register int j, i;
  Dcu = dst; Sc = src;
  Dcl = dst;
  // fill dst row and column at the same time, following src columns
  for (j = 0; j < nC; j++) {
    Dru = Dcu;
    Drl = Dcl;
    Sr = Sc;
    for (i = 0; i < j; i++) {
      *Dru = 0.0;
      *Drl = *Sr;
      Sr++;
      Dru++; 
      Drl += ldD;
    }
    // copy the diagonal entry or set to 1.0
    *Dru = unit ? 1.0 : *Sr;
    Sc += ldS;          // next column in source
    Dcu += ldD;         // next column for upper triagonal
    Dcl++;              // next row for lower triagonal
  }
}

static inline
void colcpy_trans_lower_and_zero(double *dst, int ldD, const double *src,
                                 int ldS, int nR, int nC, int unit)
{
  //assert(nR == nC);
  register double *Dcu, *Dcl, *Drl, *Dru;
  register const double *Sc, *Sr;
  register int j, i;
  Dcu = dst; Sc = src;
  Dcl = dst;
  // fill dst row and column at the same time, following src columns
  for (j = 0; j < nC; j++) {
    // start at same point and diverge down (Drl) and right (Dru)
    Dru = Dcu + j;
    Drl = Dcl + j;
    // start of data on column, j'th row (diagonal entry)
    Sr = Sc + j;
    // diagonal entry
    *Dru = unit ? 1.0 : *Sr;
    Sr++;
    Dru += ldD;
    Drl++;     
    // off diagonal entries
    for (i = 1; i < nC-j; i++) {
      *Dru = *Sr;
      *Drl = 0.0;
      Sr++;
      Dru += ldD;       // next column in row
      Drl++;            // next row in column 
    }
    // NEXT column in source
    Sc += ldS;
    // next column for upper triagonal
    Dcu += ldD;
    // next column for lower triagonal
    Dcl += ldD;
  }
}


#endif

// Local Variables:
// indent-tabs-mode: nil
// End:


