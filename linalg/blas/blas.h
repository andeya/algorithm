/*
 * These prototypes are from CVXOPT source file blas.c
 */

#ifndef BLAS_H
#define BLAS_H
#define USE_CBLAS_ZDOT 1

/* BLAS 1 prototypes */
extern void dswap_(int *n, double *x, int *incx, double *y, int *incy);
extern void zswap_(int *n, void *x, int *incx, void *y,
    int *incy);
extern void dscal_(int *n, double *alpha, double *x, int *incx);
extern void zscal_(int *n, void *alpha, void *x, int *incx);
extern void zdscal_(int *n, double *alpha, void *x, int *incx);
extern void dcopy_(int *n, double *x, int *incx, double *y, int *incy);
extern void zcopy_(int *n, void *x, int *incx, void *y,
    int *incy);
extern void daxpy_(int *n, double *alpha, double *x, int *incx,
    double *y, int *incy);
extern void zaxpy_(int *n, void *alpha, void *x, int *incx,
    void *y, int *incy);
extern double ddot_(int *n, double *x, int *incx, double *y, int *incy);
extern void zdotcsub_(int *n, void *x, int *incx, void *y, int *incy, void *result);
extern void zdotusub_(int *n, void *x, int *incx, void *y, int *incy, void *result);
extern double dnrm2_(int *n, double *x, int *incx);
extern double dznrm2_(int *n, void *x, int *incx);
extern double dasum_(int *n, double *x, int *incx);
extern double dzasum_(int *n, void *x, int *incx);

extern double dnrm2sub_(int *n, double *x, int *incx, double *result);
extern double dznrm2sub_(int *n, void *x, int *incx, void *result);
extern double dasumsub_(int *n, double *x, int *incx, double *result);
extern double dzasumsub_(int *n, void *x, int *incx, void *result);

extern int idamax_(int *n, double *x, int *incx);
extern int izamax_(int *n, void *x, int *incx);

extern int idamaxsub_(int *n, double *x, int *incx, int *result);
extern int izamaxsub_(int *n, void *x, int *incx, int *result);


/* BLAS 2 prototypes */
extern void dgemv_(char* trans, int *m, int *n, double *alpha,
    double *A, int *lda, double *x, int *incx, double *beta, double *y,
    int *incy);
extern void zgemv_(char* trans, int *m, int *n, void *alpha,
    void *A, int *lda, void *x, int *incx, void *beta,
    void *y, int *incy);
extern void dgbmv_(char* trans, int *m, int *n, int *kl, int *ku,
    double *alpha, double *A, int *lda, double *x, int *incx,
    double *beta, double *y,  int *incy);
extern void zgbmv_(char* trans, int *m, int *n, int *kl, int *ku,
    void *alpha, void *A, int *lda, void *x, int *incx,
    void *beta, void *y,  int *incy);
extern void dsymv_(char *uplo, int *n, double *alpha, double *A,
    int *lda, double *x, int *incx, double *beta, double *y, int *incy);
extern void zhemv_(char *uplo, int *n, void *alpha, void *A,
    int *lda, void *x, int *incx, void *beta, void *y,
    int *incy);
extern void dsbmv_(char *uplo, int *n, int *k, double *alpha, double *A,
    int *lda, double *x, int *incx, double *beta, double *y, int *incy);
extern void zhbmv_(char *uplo, int *n, int *k, void *alpha,
    void *A, int *lda, void *x, int *incx, void *beta,
    void *y, int *incy);
extern void dtrmv_(char *uplo, char *trans, char *diag, int *n,
    double *A, int *lda, double *x, int *incx);
extern void ztrmv_(char *uplo, char *trans, char *diag, int *n,
    void *A, int *lda, void *x, int *incx);
extern void dtbmv_(char *uplo, char *trans, char *diag, int *n, int *k,
    double *A, int *lda, double *x, int *incx);
extern void ztbmv_(char *uplo, char *trans, char *diag, int *n, int *k,
    void *A, int *lda, void *x, int *incx);
extern void dtrsv_(char *uplo, char *trans, char *diag, int *n,
    double *A, int *lda, double *x, int *incx);
extern void ztrsv_(char *uplo, char *trans, char *diag, int *n,
    void *A, int *lda, void *x, int *incx);
extern void dtbsv_(char *uplo, char *trans, char *diag, int *n, int *k,
    double *A, int *lda, double *x, int *incx);
extern void ztbsv_(char *uplo, char *trans, char *diag, int *n, int *k,
    void *A, int *lda, void *x, int *incx);
extern void dger_(int *m, int *n, double *alpha, double *x, int *incx,
    double *y, int *incy, double *A, int *lda);
extern void zgerc_(int *m, int *n, void *alpha, void *x,
    int *incx, void *y, int *incy, void *A, int *lda);
extern void zgeru_(int *m, int *n, void *alpha, void *x,
    int *incx, void *y, int *incy, void *A, int *lda);
extern void dsyr_(char *uplo, int *n, double *alpha, double *x,
    int *incx, double *A, int *lda);
extern void zher_(char *uplo, int *n, double *alpha, void *x,
    int *incx, void *A, int *lda);
extern void dsyr2_(char *uplo, int *n, double *alpha, double *x,
    int *incx, double *y, int *incy, double *A, int *lda);
extern void zher2_(char *uplo, int *n, void *alpha, void *x,
    int *incx, void *y, int *incy, void *A, int *lda);

extern void dspmv_(char *uplo, int *n, double *alpha, double *Ap, double *x,
		   int *incx, double *beta, double *y, int *incy);
extern void dspr_(char *uplo, int *n, double *alpha, double *x,
		  int *incx, double *ap);
extern void dspr2_(char *uplo, int *n, double *alpha, double *x,
		   int *incx, double *y, int *incy, double *ap);
extern void dtpmv_(char *uplo, char *transa, char *diag, int *n, double *Ap,
		   double *x, int *incx);
extern void dtpsv_(char *uplo, char *transa, char *diag, int *n, double *Ap,
		   double *x, int *incx);

/* BLAS 3 prototypes */
extern void dgemm_(char *transa, char *transb, int *m, int *n, int *k,
    double *alpha, double *A, int *lda, double *B, int *ldb,
    double *beta, double *C, int *ldc);
extern void zgemm_(char *transa, char *transb, int *m, int *n, int *k,
    void *alpha, void *A, int *lda, void *B, int *ldb,
    void *beta, void *C, int *ldc);
extern void dsymm_(char *side, char *uplo, int *m, int *n,
    double *alpha, double *A, int *lda, double *B, int *ldb,
    double *beta, double *C, int *ldc);
extern void zsymm_(char *side, char *uplo, int *m, int *n,
    void *alpha, void *A, int *lda, void *B, int *ldb,
    void *beta, void *C, int *ldc);
extern void zhemm_(char *side, char *uplo, int *m, int *n,
    void *alpha, void *A, int *lda, void *B, int *ldb,
    void *beta, void *C, int *ldc);
extern void dsyrk_(char *uplo, char *trans, int *n, int *k,
    double *alpha, double *A, int *lda, double *beta, double *B,
    int *ldb);
extern void zsyrk_(char *uplo, char *trans, int *n, int *k,
    void *alpha, void *A, int *lda, void *beta, void *B,
    int *ldb);
extern void zherk_(char *uplo, char *trans, int *n, int *k,
    double *alpha, void *A, int *lda, double *beta, void *B,
    int *ldb);
extern void dsyr2k_(char *uplo, char *trans, int *n, int *k,
    double *alpha, double *A, int *lda, double *B, int *ldb,
    double *beta, double *C, int *ldc);
extern void zsyr2k_(char *uplo, char *trans, int *n, int *k,
    void *alpha, void *A, int *lda, void *B, int *ldb,
    void *beta, void *C, int *ldc);
extern void zher2k_(char *uplo, char *trans, int *n, int *k,
    void *alpha, void *A, int *lda, void *B, int *ldb,
    double *beta, void *C, int *ldc);
extern void dtrmm_(char *side, char *uplo, char *transa, char *diag,
    int *m, int *n, double *alpha, double *A, int *lda, double *B,
    int *ldb);
extern void ztrmm_(char *side, char *uplo, char *transa, char *diag,
    int *m, int *n, void *alpha, void *A, int *lda, void *B,
    int *ldb);
extern void dtrsm_(char *side, char *uplo, char *transa, char *diag,
    int *m, int *n, double *alpha, double *A, int *lda, double *B,
    int *ldb);
extern void ztrsm_(char *side, char *uplo, char *transa, char *diag,
    int *m, int *n, void *alpha, void *A, int *lda, void *B,
    int *ldb);

#endif
