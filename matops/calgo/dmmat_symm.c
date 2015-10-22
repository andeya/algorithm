
// Copyright (c) Harri Rautila, 2012,2013

// This file is part of github.com/hrautila/matops package. It is free software,
// distributed under the terms of GNU Lesser General Public License Version 3, or
// any later version. See the COPYING tile included in this archive.

#include <stdio.h>

#include "cmops.h"
#include "colcpy.h"

extern void
_dblock_ddot_sse(double *Cc, const double *Aroot, const double *Bc, double alpha,
                 int ldC, int ldA, int ldB, int nSL, int nRE, int nVP);


void _dblock_symm_cpy(mdata_t *C, const mdata_t *A, const mdata_t *B,
                      double alpha, double beta, int flags,
                      int nP, int S, int L, int R, int E, int vlen)
{
  int j, k, vpS, vpL, nC, nB, nA;
  const double *Bc, *Ac, *AvpS;
  //const double *Br0, *Br1, *Br2, *Br3;
  double *Cc; //, *c0, *c1, *c2, *c3;
  //double Cpy[MAX_NB_DDOT*MAX_MB_DDOT]  __attribute__((aligned(16)));
  double Acpy[MAX_VP_DDOT*MAX_MB_DDOT] __attribute__((aligned(16)));
  double Bcpy[MAX_VP_DDOT*MAX_NB_DDOT] __attribute__((aligned(16)));
  int unit = flags & MTX_UNIT ? 1 : 0;


  if (vlen > nP || vlen <= 0) {
    vlen = nP;
  }

  //printf("0. nP=%d, L=%d, S=%d, E=%d, R=%d, vlen=%d\n", nP, L, S, E, R, vlen);

  // row stride in Cpy 
  //nC = E - R;
  //nC += (nC & 0x1); // increments by 1 if not even.

  // Copy C block to local buffer
  Cc = &C->md[S*C->step+R];
  //colcpy(Cpy, nC, Cc, C->step, E-R, L-S);
  nC = C->step;

  // TODO: scaling with beta ....
  dscale_tile(Cc, nC, beta, E-R, L-S);

  // nP is columns in A, 
  vpS = 0;
  vpL = vlen < R ? vlen : R;

  // 1. this is the panel left of diagonal block (LOWER) or above the diagonal (UPPER) 
  //    vps, vpL < R: work it like notrans
  while (vpS < R) {
    nB = vpL-vpS;
    nB += (nB & 0x1);
    nA = nB;
    //printf("1. vpS=%d, vpL=%d, nC=%d, nB=%d, nA=%d\n", vpS, vpL, nC, nB, nA);

    // viewport starts in B, A
    if (flags & MTX_LOWER) {
      // The panel left of diagonal block
      Bc = &B->md[S*B->step + vpS];
      AvpS = &A->md[vpS*A->step + R];
      // transpose A on copy to be able to DOT operations.
      colcpy(Bcpy, nB, Bc, B->step, vpL-vpS, L-S);
      colcpy4_trans(Acpy, nA, AvpS, A->step, E-R, vpL-vpS);
    } else {
      // The panel above of diagonal block
      Bc = &B->md[S*B->step + vpS];
      AvpS = &A->md[R*A->step + vpS];

      colcpy(Bcpy, nB, Bc, B->step, vpL-vpS, L-S);
      colcpy(Acpy, nA, AvpS, A->step, vpL-vpS, E-R); //, vpL-vpS);
    }
    //printf("1. update: A=\n"); print_tile(Acpy, nA, vpL-vpS, E-R);
    //printf("1. update: B=\n"); print_tile(Bcpy, nB, vpL-vpS, L-S);

    if (flags & MTX_LEFT) {
      _dblock_ddot_sse(Cc, Acpy, Bcpy, alpha, nC, nA, nB, L-S, E-R, vpL-vpS);
    } else {
    }
    //printf("1. post update: C=\n"); print_tile(Cpy, nC, E-R, L-S);

    vpS = vpL;
    vpL += vlen;
    if (vpL > R) {
      vpL = R;
    }
  }
  
  // 2. this is the diagonal block, with upper part untouchable
  //    R <= vps, vpL < E: diagonal part, copy_and_fill
  //    here vpS == R, update vpL += E - R as this block is square, diagonal
  vpL += E-R;
  while (vpS < E) {
    nB = vpL-vpS;
    nB += (nB & 0x1);
    nA = nB;

    //printf("2. vpS=%d, vpL=%d, nC=%d, nB=%d, nA=%d\n", vpS, vpL, nC, nB, nA);
    // viewport starts in B, A
    if (flags & MTX_LOWER) {
      // upper part of source untouchable, copy diagonal block and fill upper part
      Bc = &B->md[S*B->step + vpS];
      AvpS = &A->md[vpS*A->step + R];
      //print_tile(Acpy, nA, E-R, E-R);
      colcpy(Bcpy, nB, Bc, B->step, vpL-vpS, L-S);
      colcpy_fill_up(Acpy, nA, AvpS, A->step, E-R, E-R, unit);
    } else {
      // lower part of source untouchable, copy diagonal block and fill lower part
      Bc = &B->md[S*B->step + vpS];
      AvpS = &A->md[vpS*A->step + R];
      colcpy(Bcpy, nB, Bc, B->step, vpL-vpS, L-S);
      colcpy_fill_low(Acpy, nA, AvpS, A->step, E-R, E-R, unit);
    }

    if (flags & MTX_LEFT) {
      _dblock_ddot_sse(Cc, Acpy, Bcpy, alpha, nC, nA, nB, L-S, E-R, vpL-vpS);
    }
    //printf("2. post update: C=\n"); print_tile(Cpy, nC, E-R, L-S);

    vpS = vpL;
    vpL += vlen;
    if (vpL > E) {
      vpL = E;
    }
  }

  // 3. this is rest of the panel rows below or right of diagonal block.
  //    vps, vpL >= E && < nP: rest of the row like transA case.
  //    here vpS == E, and vpL == vpS + vlen
  vpL = vpS + vlen;
  if (vpL > nP) {
    vpL = nP;
  }
  while (vpS < nP) {
    nB = vpL-vpS;
    nB += (nB & 0x1);
    nA = nB;

    //printf("3. vpS=%d, vpL=%d, nC=%d, nB=%d, nA=%d\n", vpS, vpL, nC, nB, nA);
    if (flags & MTX_LOWER) {
      // this is rest of the panel rows below of the diagonal block.
      Bc = &B->md[S*B->step + vpS];
      AvpS = &A->md[R*A->step + vpS];
      colcpy(Bcpy, nB, Bc, B->step, vpL-vpS, L-S);
      colcpy(Acpy, nA, AvpS, A->step, vpL-vpS, E-R);
    } else {
      // this is rest of the panel rows right of the diagonal block.
      Bc = &B->md[S*B->step + vpS];
      AvpS = &A->md[vpS*A->step + R];

      // transpose A on copy to be able to DOT operations.
      colcpy(Bcpy, nB, Bc, B->step, vpL-vpS, L-S);
      colcpy4_trans(Acpy, nA, AvpS, A->step, E-R, vpL-vpS);
    }

    if (flags & MTX_LEFT) {
      // C += alpha * A * B
      _dblock_ddot_sse(Cc, Acpy, Bcpy, alpha, nC, nA, nB, L-S, E-R, vpL-vpS);
    } else {
      // C += alpha * B * A
    }
    //printf("3. post update: C=\n"); print_tile(Cpy, nC, E-R, L-S);

    vpS = vpL;
    vpL += vlen;
    if (vpL > nP) {
      vpL = nP;
    }
  }
  // copy back.
  //colcpy(Cc, C->step, Cpy, nC, E-R, L-S);
}


void dmult_symm_blocked(mdata_t *C, const mdata_t *A, const mdata_t *B,
                        double alpha, double beta, int flags,
                        int P, int S, int L, int R, int E,
                        int vlen, int NB, int MB)
{
  int i, j, nI, nJ;

  if (L-S <= 0 || E-R <= 0) {
    return;
  }

  // restrict block sizes as data is copied to aligned buffers of predefined max sizes.
  if (NB > MAX_NB_DDOT || NB <= 0) {
    NB = MAX_NB_DDOT;
  }
  if (MB > MAX_MB_DDOT || MB <= 0) {
    MB = MAX_MB_DDOT;
  }
  if (vlen> MAX_VP_DDOT || vlen <= 0) {
    vlen = MAX_VP_DDOT;
  }

  // A is square matrix: 
  for (j = S; j < L; j += NB) {
    nJ = L - j < NB ? L - j : NB;
    for (i = R; i < E; i += MB) {
      nI = E - i < MB ? E - i : MB;
      _dblock_symm_cpy(C, A, B, alpha, beta, flags, P, j, j+nJ, i, i+nI, vlen);
      //printf("\nC=\n"); print_tile(C->md, C->step, E-R, L-S);
    }
  }
}



// C += A*B; A is the diagonal block
void _dblock_mult_diag(mdata_t *C, const mdata_t *A, const mdata_t *B,
                       double alpha, int flags, 
                       int nP, int nSL, int nRE, int vlen, cbuf_t *Acpy, cbuf_t *Bcpy)
{
  // assert (nSL == nRE)
  int unit = flags & MTX_UNIT ? 1 : 0;
  int nA, nB, nAC;

  if (nP == 0)
    return;

  nA = nRE + (nRE & 0x1);
  nB = nA;
  nAC = flags & MTX_RIGHT ? nSL : nRE;
  
  if (flags & MTX_LOWER) {
    // upper part of source untouchable, copy diagonal block and fill upper part
    colcpy_fill_up(Acpy->data, nA, A->md, A->step, nAC, nAC, unit);
  } else {
    // lower part of source untouchable, copy diagonal block and fill lower part
    colcpy_fill_low(Acpy->data, nA, A->md, A->step, nAC, nAC, unit);
  }

  if (flags & MTX_RIGHT) {
    colcpy_trans(Bcpy->data, nB, B->md, B->step, nRE, nSL);
  } else {
    colcpy(Bcpy->data, nB, B->md, B->step, nRE, nSL);
  }
  if (flags & MTX_RIGHT) {
    _dblock_ddot_sse(C->md, Bcpy->data, Acpy->data, alpha, C->step, nB, nA, nSL, nRE, nP);
  } else {
    _dblock_ddot_sse(C->md, Acpy->data, Bcpy->data, alpha, C->step, nA, nB, nSL, nRE, nP);
  }
  //printf("2. post update: C=\n"); print_tile(Cpy, nC, E-R, L-S);

}

void dmult_symm_blocked2(mdata_t *C, const mdata_t *A, const mdata_t *B,
                         double alpha, double beta, int flags,
                         int P, int S, int L, int R, int E,
                         int vlen, int NB, int MB)
{
  int i, j, nI, nJ, flags1, flags2;
  mdata_t A0, B0, C0, *Ap, *Bp;
  double Abuf[MAX_VP_ROWS*MAX_VP_COLS] __attribute__((aligned(64)));
  double Bbuf[MAX_VP_ROWS*MAX_VP_COLS] __attribute__((aligned(64)));
  cbuf_t Acpy = {Abuf, MAX_VP_ROWS*MAX_VP_COLS};
  cbuf_t Bcpy = {Bbuf, MAX_VP_ROWS*MAX_VP_COLS};

  if (L-S <= 0 || E-R <= 0) {
    return;
  }

  // restrict block sizes as data is copied to aligned buffers of predefined max sizes.
  if (NB > MAX_VP_COLS || NB <= 0) {
    NB = MAX_VP_COLS;
  }
  if (MB > MAX_VP_ROWS || MB <= 0) {
    MB = MAX_VP_ROWS;
  }
  if (vlen > MAX_VP_ROWS || vlen <= 0) {
    vlen = MAX_VP_ROWS;
  }

  C0.step = C->step;
  A0.step = A->step;
  B0.step = B->step;

  flags1 = 0;
  flags2 = 0;

  if (flags & MTX_LEFT) {
    /*
      P is A, B common dimension, e.g. P cols in A and P rows in B.
    
      [R,R] [E,E] define block on A diagonal that divides A in three blocks
      if A is upper:
        A0 [0, R] [R, E]; B0 [0, S] [R, L] (R rows,cols in P); (A transposed)
        A1 [R, R] [E, E]; B1 [R, S] [E, L] (E-R rows,cols in P)
        A2 [R, E] [E, N]; B2 [E, S] [N, L] (N-E rows, cols in  P)
      if A is LOWER:
        A0 [R, 0] [E, R]; B0 [0, S] [R, L]
        A1 [R, R] [E, E]; B1 [R, S] [E, L] (diagonal block, fill_up);
        A2 [E, R] [E, N]; B2 [E, S] [N, L] (A transpose)
        
      C = A0*B0 + A1*B1 + A2*B2
    */
    flags1 |= flags & MTX_UPPER ? MTX_TRANSA : 0;
    flags2 |= flags & MTX_LOWER ? MTX_TRANSA : 0;

    for (i = R; i < E; i += MB) {
      nI = E - i < MB ? E - i : MB;

      // for all column of C, B ...
      for (j = S; j < L; j += NB) {
        nJ = L - j < NB ? L - j : NB;
        C0.md = &C->md[j*C->step + i];
      
        dscale_tile(C0.md, C0.step, beta, nI, nJ);

        // above|left diagonal
        A0.md = flags & MTX_UPPER ? &A->md[i*A->step] : &A->md[i];
        B0.md = &B->md[j*B->step];
        //_dblock_mult_panel(&C0, &A0, &B0, alpha, flags1, i, nJ, nI, vlen, &Acpy, &Bcpy);
        _dmult_mm_intern(&C0, &A0, &B0, alpha, flags1, i, nJ, nI, vlen, NB, MB, &Acpy, &Bcpy);

        // diagonal block
        A0.md = &A->md[i*A->step + i];
        B0.md = &B->md[j*B->step + i];
        _dblock_mult_diag(&C0, &A0, &B0, alpha, flags, nI, nJ, nI, vlen, &Acpy, &Bcpy);

        // right|below of diagonal
        A0.md = flags & MTX_UPPER ? &A->md[(i+nI)*A->step + i] : &A->md[i*A->step + i+nI];
        B0.md = &B->md[j*B->step + i+nI];
        //_dblock_mult_panel(&C0, &A0, &B0, alpha, flags2, P-i-nI, nJ, nI, vlen, &Acpy, &Bcpy);
        _dmult_mm_intern(&C0, &A0, &B0, alpha, flags2, P-i-nI, nJ, nI, vlen, NB, MB, &Acpy, &Bcpy);
      }
    }
  } else {

    /*
      P is A, B common dimension, e.g. P cols in A and P rows in B.
    
      C = B * A;
      [S,S] [L,L] define block on A diagonal that divides A in three blocks
      if A is upper:
        A0 [0, S] [S, S]; B0 [R, 0] [E, S] (R rows,cols in P); (A transposed)
        A1 [S, S] [L, L]; B1 [R, S] [E, L] (E-R rows,cols in P)
        A2 [S, L] [L, N]; B2 [R, L] [E, N] (N-E rows, cols in  P)
      if A is LOWER:
        A0 [S, 0] [S, S]; B0 [R, 0] [E, S]
        A1 [S, S] [L, L]; B1 [R, S] [E, L] (diagonal block, fill_up);
        A2 [L, S] [N, L]; B2 [R, L] [E, N] (A transpose)
        
      C = A0*B0 + A1*B1 + A2*B2
    */

    register int nR, nC, ic, ir;
    flags1 = flags & MTX_TRANSB ? MTX_TRANSA : 0;
    flags2 = flags & MTX_TRANSB ? MTX_TRANSA : 0;

    flags1 |= flags & MTX_LOWER ? MTX_TRANSB : 0;
    flags2 |= flags & MTX_UPPER ? MTX_TRANSB : 0;

    for (ic = S; ic < L; ic += NB) {
      nC = L - ic < NB ? L - ic : NB;

      // for all rows of C, B ...
      for (ir = R; ir < E; ir += MB) {
        nR = E - ir < MB ? E - ir : MB;

        C0.md = &C->md[ic*C->step + ir];
      
        dscale_tile(C0.md, C0.step, beta, nR, nC);

        // above|left diagonal
        A0.md = flags & MTX_UPPER ? &A->md[ic*A->step] : &A->md[ic];
        B0.md = &B->md[ir];
        //_dblock_mult_panel(&C0, &B0, &A0, alpha, flags1, ic, nC, nR, vlen, &Acpy, &Bcpy);
        _dmult_mm_intern(&C0, &B0, &A0, alpha, flags1, ic, nC, nR, vlen, NB, MB, &Acpy, &Bcpy);

        // diagonal block
        A0.md = &A->md[ic*A->step + ic];
        B0.md = &B->md[ic*B->step+ir];
        _dblock_mult_diag(&C0, &A0, &B0, alpha, flags, nC, nC, nR, vlen, &Acpy, &Bcpy);

        // right|below of diagonal
        A0.md = flags & MTX_UPPER ? &A->md[(ic+nC)*A->step + ic] : &A->md[ic*A->step +ic+nC];
        B0.md = &B->md[(ic+nC)*B->step+ir];
        //_dblock_mult_panel(&C0, &B0, &A0, alpha, flags2, P-ic-nC, nC, nR, vlen, &Acpy, &Bcpy);
        _dmult_mm_intern(&C0, &B0, &A0, alpha, flags2, P-ic-nC, nC, nR, vlen, NB, MB, &Acpy, &Bcpy);
      }
    }
  }
}




// Local Variables:
// indent-tabs-mode: nil
// End:
