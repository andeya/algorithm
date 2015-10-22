/*
 * interfaceBLAS-LAPACK.h
 *
 *  Created on: Jun 22, 2010
 *      Author: awblocker
 */

#ifndef INTERFACEBLASLAPACK_H_
#define INTERFACEBLASLAPACK_H_

void dscal(int N, double ALPHA, double * X, int INCX);

double dasum(int N, double * X, int INCX);

void dcopy(int N, double * X, int INCX, double * Y, int INCY);

void daxpy(int N, double ALPHA, double * X, int INCX, double * Y, int INCY);

void dgemv(char TRANS, int M, int N, double ALPHA, double * A, int LDA,
    double * X, int INCX, double BETA, double * Y, int INCY);

void dsyrk(char UPLO, char TRANS, int N, int K, double ALPHA,
    double* A, int LDA, double BETA, double* C, int LDC);

int dposv(char UPLO, int N, int NRHS, double * A, int LDA,
    double * B, int LDB);

int dgels(char TRANS, int M, int N, int NHRS, double * A, int LDA,
    double * B, int LDB);

#endif /* INTERFACEBLASLAPACK_H_ */
