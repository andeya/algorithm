
// Copyright (c) Harri Rautila, 2012

// This file is part of github.com/hrautila/gblas package. It is free software,
// distributed under the terms of GNU Lesser General Public License Version 3, or
// any later version. See the COPYING tile included in this archive.


#define VLEN_DEFAULT 30

// Matrix data is assumed to be in column major order. 
// Matrix C block defined by rows R:E, columns S:L,
// matrix B column panel by S, L; matrix  A rows panel by R:E.
void matmult_block_notrans(double *C, const double *A, const double *B,
                           const double alpha,
                           int M, int N, int P, int S, int L, int R, int E)
{
  // M is rows in C, A
  // P is cols in A, rows in B
  int i, j, k;
  const double *Br, *Ac, *Bc;
  double *Cc;
  register const double *Ar;
  register double *Cr;
  register double coeff;

  Cc = &C[S*M+R];	// block start C[R, S]
  Bc = &B[S*P];		// column panel start B[:,S]
  Ac = &A[R];		// row panel start A[R,:]
  for (j = S; j < L; j++) {
    Ac = A + R;
    Br = Bc;
    for (k = 0; k < P; k++) {
      Cr = Cc; 
      Ar = Ac;
       if (*Br != 0.0) {
	// C[:,j] += A[:,k]*B[k,j] 
	coeff = *Br * alpha; 
	for (i = R; i < E; i++) {
	  *Cr += (*Ar) * coeff;
	  Cr++;
          Ar++; 
	}
      }
      // move to next row in B, here ar points to first row of next column in A
      Br++; 
      Ac += M; // next column in A 
    }
    // forward to start of next column in C 
    Cc += M;
    Bc += P;
  }
}


// C is M*N, A is M*P and B is P*N column major ordered matrix data elements.
// S:L are column indexes for B column panel and C block. R:E are row indexes
// for A row panel and C block. Parameter vlen is 'viewport' length A columns
// and B rows. A and B panel data is accumulated to C in 'vlen' blocks in order
// not to access A, B panels for full length (P) and receiving better cache
// hit ratio for A elements that are accessed multiple times. If R:E*vlen is small
// enought it maybe fits into L1 cache. (Not very scientific!)
// 
// Index value ranges:
//      0 <= S < L < N
//      0 <= R < E < M
//      0 < vlen < P
//
void matmult_vp_notrans(double *C, const double *A, const double *B, double alpha,
                        int M, int N, int P, int S, int L, int R, int E, int vlen)
{
  int j, k, vpS, vpL;
  const double *Bc, *Br, *Ac, *AvpS;
  double *Cc;
  register const double *Ar;
  register int i;
  register double *Cr;
  register double coeff;

  if (vlen == 0) {
    vlen = VLEN_DEFAULT;
  }
  vpS = 0;
  vpL = vlen < P ? vlen : P;
  while (vpS < P) {
    Cc = &C[S*M+R];             // block start C[R, S]
    Bc = &B[S*P + vpS];         // column viewport start in panel B
    AvpS = &A[vpS*M + R];       // row viewport start in panel A

    for (j = S; j < L; j++) {
      Ac = AvpS;      // reset A to start of viewport
      Br = Bc;
      for (k = vpS; k < vpL; k++) {
        // move A,C indexes to first row in current column 
        Cr = Cc;
        Ar = Ac;
        if (*Br != 0.0) {
          coeff = (*Br) * alpha;
          // TODO: some loop unrolling here???
          for (i = R; i < E; i++) {
            *Cr += (*Ar) * coeff;
            Cr++;
            Ar++;
          }
        } 
        // move to next row in B, next column in A
        Br++;
        Ac += M;
      }
      // forward to start of next column in C, B
      Cc += M;
      Bc += P;
    }
    
    vpS = vpL;
    vpL += vlen;
    if (vpL > P) {
      vpL = P;
    }
  }
}


// Same as above with inner most loop partially unrolled. Speed up ~8%.
void matmult_vpur_notrans(double *C, const double *A, const double *B, double alpha,
                          int ldC, int ldA, int ldB,
                          int M, int N, int P, int S, int L, int R, int E, int vlen)
{
  int j, k, vpS, vpL;
  int SLcnt, REcnt, REmod;
  const double *Bc, *Br, *Ac, *AvpS;
  double *Cc;
  register const double *Ar;
  register int i;
  register double *Cr;
  register double coeff, f0, f1, f2, f3;

  SLcnt = L - S;
  REcnt = E - R;
  REmod = REcnt % 4;
  REcnt -= REmod;
  if (vlen == 0) {
    vlen = 32;
  }
  vpS = 0;
  vpL = vlen < P ? vlen : P;
  while (vpS < P) {
    // block start C[R, S]
    Cc = &C[S*ldC+R];
    // column viewport start in panel B[:,S]
    Bc = &B[S*ldB + vpS];
    AvpS = &A[vpS*ldA + R];

    for (j = S; j < L; j++) {
      Ac = AvpS;       // row viewport start A[R,:]
      Br = Bc;
      for (k = vpS; k < vpL; k++) {
        Cr = Cc;        // move C index to first row in current column 
        Ar = Ac;
        coeff = *Br;
        if (coeff != 0.0) {
          coeff *= alpha;
          for (i = 0; i < REmod; i++) {
            f0 = Ar[0] * coeff;
            Cr[0] += f0;
            Cr++;
            Ar++;
          }
          for (i = 0; i < REcnt; i += 4) {
            f0 = Ar[0] * coeff;
            Cr[0] += f0;
            f1 = Ar[1] * coeff;
            Cr[1] += f1;
            f2 = Ar[2] * coeff;
            Cr[2] += f2;
            f3 = Ar[3] * coeff;
            Cr[3] += f3;
            Cr += 4;
            Ar += 4;
          }
        } 
        // move to next row in B, next column in A
        Br++;
        Ac += ldA;
      }
      // forward to start of next column in C, B
      Cc += ldC;
      Bc += ldB;
    }
    
    vpS = vpL;
    vpL += vlen;
    if (vpL > P) {
      vpL = P;
    }
  }
  return;
}



// Local Variables:
// indent-tabs-mode: nil
// End:
