/*
 * wls.c
 *
 *  Created on: Jun 18, 2010
 *      Author: awblocker
 */

#include <math.h>
#include "interfaceBLAS-LAPACK.h"

int wls(double* X, int n, int p, double* y, double* w, char method,
        double* XTX, double *sqw, double* sqwX, double* sqwy, double* coef) {
    // Assuming column-major order
    // Initializations
    int i, j;

    // Compute sqw
    for (i=0; i<n; i++)
    {
        sqw[i] = sqrt(w[i]);

        // Compute sqwX and sqwy
        for (j=0; j<p; j++) {
            sqwX[i + j*n] = X[i + j*n] * sqw[i];
        }
        sqwy[i] = sqw[i] * y[i];
    }

    // Method selection; q -> QR, c (or anything else) -> Cholesky
    int info;
    if (method == 'q') {
        // Obtain least-squares coefficients via LAPACK DGELS driver routine
        info = dgels('n', n, p, 1, sqwX, n, sqwy, n);

        // Copy solutions to coef vector
        dcopy(p, sqwy, 1, coef, 1);
    } else {
        // Compute XTX
        dsyrk('u', 't', p, n, 1, sqwX, n, 0, XTX, p);

        // Compute Xy
        dgemv('t', n, p, 1, sqwX, n, sqwy, 1, 0, coef, 1);

        // Obtain least-squares coefficients by solving normal equations
        info = dposv('u', p, 1, XTX, p, coef, p);
    }

    return info;
}

int calcFitted(double* X, int n, int p,
        double* y,
        double* coef,
        double* fitted)
{
    // X \beta = \hat{X}
    dgemv('n', n, p, 1, X, n, coef, 1, 0, fitted, 1);

    return 0;
}

int calcResid(double* X, int n, int p,
        double* y,
        double* coef,
        double* resid)
{
    // Calculate fitted values; store temporarily in resid
    dgemv('n', n, p, 1, X, n, coef, 1, 0, resid, 1);

    // Calculate residuals
    dscal(n, -1, resid, 1);
    daxpy(n, 1, y, 1, resid, 1);

    return 0;
}
