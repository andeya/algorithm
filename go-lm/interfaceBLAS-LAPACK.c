/*
 * interfaceBLAS-LAPACK.c
 *
 *  Created on: Jun 22, 2010
 *      Author: awblocker
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

inline int max(int a, int b) {
  return a > b ? a : b;
}

extern void dscal_(int * N, double * ALPHA, double * X, int * INCX);

void dscal(int N, double ALPHA, double * X, int INCX)
{
  dscal_(&N, &ALPHA, X, &INCX);
}

extern double dasum_(int * N, double * X, int * INCX);

double dasum(int N, double * X, int INCX)
{
  return dasum_(&N, X, &INCX);
}

extern void dcopy_(int * N, double * X, int * INCX, double * Y, int * INCY);

void dcopy(int N, double * X, int INCX, double * Y, int INCY)
{
  dcopy_(&N, X, &INCX, Y, &INCY);
}

extern void daxpy_(int * N, double * ALPHA, double * X, int * INCX,
    double * Y, int * INCY);

void daxpy(int N, double ALPHA, double * X, int INCX, double * Y, int INCY)
{
  daxpy_(&N, &ALPHA, X, &INCX, Y, &INCY);
}

extern void dgemv_(char * TRANS, int * M, int * N, double * ALPHA, double * A,
    int * LDA, double * X, int * INCX, double * BETA,
    double * Y, int * INCY);

void dgemv(char TRANS, int M, int N, double ALPHA, double * A, int LDA,
    double * X, int INCX, double BETA, double * Y, int INCY)
{
  dgemv_(&TRANS, &M, &N, &ALPHA, A, &LDA, X, &INCX, &BETA, Y, &INCY);
}

extern void dsyrk_(char * UPLO, char * TRANS, int * N, int * K,
    double * ALPHA, double * A, int * LDA, double * BETA,
    double * C,  int *LDC);

void dsyrk(char UPLO, char TRANS, int N, int K, double ALPHA,
    double* A, int LDA, double BETA, double* C, int LDC)
{
  dsyrk_(&UPLO, &TRANS, &N, &K, &ALPHA, A, &LDA, &BETA, C, &LDC);
}

extern void dposv_(char * UPLO, int * N, int * NRHS, double * A, int * LDA,
    double * B, int * LDB, int * INFO);

int dposv(char UPLO, int N, int NRHS, double * A, int LDA,
    double * B, int LDB)
{
  int INFO;
  dposv_(&UPLO, &N, &NRHS, A, &LDA, B, &LDB, &INFO);
  return INFO;
}

extern void dgels_(char * TRANS, int * M, int * N, int * NHRS,
    double * A, int * LDA, double * B, int * LDB, double * WORK,
    int * LWORK, int * INFO);

int dgels(char TRANS, int M, int N, int NRHS, double * A, int LDA,
    double * B, int LDB)
{
  // Set block size
  const int NB=1024;

  // Allocate work arrays of appropriate dimensions
  int LWORK;
  LWORK = M*N + max(M*N, NRHS)*NB;

  double * WORK;
  WORK = malloc(LWORK * sizeof(double));
  if (WORK==NULL)
  {
    fprintf(stderr, "Error -- out of memory\n");
    exit(1);
  }

  // Calculate least-squares solution
  int INFO;
  dgels_(&TRANS, &M, &N, &NRHS, A, &LDA, B, &LDB, WORK, &LWORK, &INFO);

  return INFO;
}

