matops
======

Matrix operations. Almost complete implementation of BLAS level 1, 2 and 3 routines for double precision floating point. All computation is in place. The implementation supports matrix views (submatrices of larger matrices) and parallel execution of matrix operations in multiple threads. Some functions already allow parallel execution for better performance.

**NOTE:** I am not working on this code base anymore. I have moved most of the functionality to a new package
**github.com/hrautila/gomas**. See there for more details.

Supported functionality is:

  Blas level 3

    Mult(C, A, B, alpha, beta, flags)           General matrix-matrix multiplication  (GEMM)
    MultSymm(C, A, B, alpha, beta, flags)       Symmetric matrix-matrix multipication (SYMM)
    MultTrm(B, A, alpha, flags)                 Triangular matrix-matrix multiplication (TRMM)  
    SolveTrm(B, A, alpha, flags)                Triangular solve with multiple RHS (TRSM)
    RankUpdateSym(C, A, alpha, beta,flags)      Symmetric matrix rank-k update (SYRK)
    RankUpdate2Sym(C, A, B, alpha, beta, flags) Symmetric matrix rank-2k update (SYR2K)

  Blas level 2

    MVMult(X, A, Y, alpha, beta, flags)         General matrix-vector multiplication (GEMV)
    MVRankUpdate(A, X, Y, alpha, flags)         General matrix rank update (GER)
    MVRankUpdateSym(A, X, alpha, flags)         Symmetric matrix rank update (SYR)
    MVRankUpdate2Sym(A, X, Y, alpha, flags)     Symmetric matrix rank 2 update (SYR2)
    MVSolveTrm(X, A, alpha, flags)              Triangular solve (TRSV)
    MVMultTrm(X, A, flags)                      Triangular matrix-vector multiplication (TRMV)

  Blas level 1

    Norm2(X, Y)         Vector norm (NRM2)
    Dot(X, Y)           Inner product (DOT)
    Swap(X, Y)          Vector-vector swap (SWAP)
    InvScale(X, alpha)  Inverse scaling of X 
    Scale(X, alpha)     Scaling of X (SCAL)

  Additional

    ScalePlus(A, B, alpha, beta, flags)         Calculate A = alpha*op(A) + beta*op(B)
    NormP(X, norm)                              Matrix or vector norm, _1, _2, _Inf
    UpdateTrm(C, A, B, alpha, beta, flags)      Triangular matrix update
    MVUpdateTrm(C, X, Y, alpha, flags)          Triangular matrix update with vectors.

  Lapack
  
    DecomposeCHOL(A, nb)                Cholesky factorization (DPOTRF)
    DecomposeLDLnoPiv(A, nb)            LDL factorization without pivoting
    DecomposeLDL(A, W, ipiv, flgs, nb)  LDL factorization with pivoting
    DecomposeLUnoPiv(A, nb)             LU factorization without pivoting
    DecomposeLU(A, pivots, nb)          LU factorization with pivoting (DGETRF)
    DecomposeQR(A, tau, nb)             QR factorization (DGEQRF)
    DecomposeQRT(A, T, W, nb)           QR factorization, compact WY version (DGEQRT)
    MultQ(C, A, tau, W, flgs, nb)       Multiply by Q  (DORMQR)
    MultQT(C, A, T, W, flgs, nb)        Multiply by Q, compact WY version (DGEMQRT)
    BuildQ(A, tau, W, nb)               Build matrix Q with ortonormal columns (DORGQR)
    BuildQT(A, T, W, nb)                Build matrix Q with ortonormal columns 
    BuildT(T, A, tau)                   Build block reflector T from elementary reflectors (DLARFT)
    SolveCHOL(B, A, flags)              Solve Cholesky factorized linear system (DPOTRS)
    SolveLDL(B, A, pivots, flags)       Solve LDL factorized linear system
    SolveLU(B, A, pivots, flags)        Solve LU factorized linear system (DGETRS)
    SolveQR(B, A, tau, W, flgs, nb)     Solve least square problem when m >= n (DGELS)
    SolveQRT(B, A, T, W, flgs, nb)      Solve least square problem when m >= n, compact WY (DGELS)
    InverseTrm(A, flags, nb)            Inverse triangular matrix (DTRTRI)

  Support functions

    TriL(A)                   Make A triangular, lower 
    TriLU(A)                  Make A triangular, lower, unit-diagonal 
    TriU(A)                   Make A triangular, upper 
    TriUU(A)                  Make A triangular, upper, unit-diagonal 

  Parameter functions

    BlockingParams(m,n,k)     Blocking size parameters for low-level functions
    NumWorkers(nwrk)          Number of threads for use in operations
    DecomposeBlockSize(nb)    Block size for blocked decomposition algorithms

This is still WORK IN PROGRESS. Consider this as beta level code, at best. 

Overall performance is compareable to ATLAS BLAS library. Some performance testing programs are in test subdirectory. Running package and performace tests requires github.com/hrautila/linalg packages as results are compared to existing BLAS/LAPACK implementation.

See the Wiki pages for some additional information. 
