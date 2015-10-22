
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "interfaceBLAS-LAPACK.h"
#include "wls.h"

// Compute sum of log-density for normal
// Omits normalizing constants
double dnorm_log(double * x, int n, double location, double scale)
{
    int i;
    double logDensity = 0, z, logScale;

    logScale = log(scale);

    for (i=0; i<n; i++)
    {
        z = (x[i] - location) / scale;
        logDensity += -0.5*z*z - logScale;
    }

    return logDensity;
}

// Compute sum of log-density for t distribution
// Omits normalizing constants
double dt_log(double * x, int n, double df, double location, double scale)
{
    int i;
    double logDensity = 0, z, logScale;

    logScale = log(scale);

    for (i=0; i<n; i++)
    {
        z = (x[i] - location) / scale;
        logDensity += -(df+1)/2 * log(1 + z*z/df) - logScale;
    }

    return logDensity;
}

/*
 * Function to run PX-EM algorithm for linear regression with t_nu residuals
 * See Meng & van Dyk 1997, JRSS-B for details
 * Returns number of iterations run
 * Log-posterior, log-likelihood, and coefficients are returned by reference
 */

int lmT(double * X, int n, int p,
        double * y,
        double nu,
        int maxIter, double tol, char method,
        double * logLikelihood,
        double * coef,
        double * tau)
{
    // Initialize workspace variables
    int iter, i;
    double logLikelihood_tm1, delta;

    /*
     * Allocate workspace arrays
     */

    // Allocate workspace matrices
    double * sqwX, * XTX;
    sqwX = malloc(n * p * sizeof(double));
    if (sqwX==NULL)
    {
        fprintf(stderr, "Error -- out of memory\n");
        exit(1);
    }

    XTX = malloc(p * p * sizeof(double));
    if (XTX==NULL)
    {
        fprintf(stderr, "Error -- out of memory\n");
        exit(1);
    }

    // Allocate workspace vectors
    double * w, * sqw, * sqwy, * resid;

    w = calloc(n, sizeof(double));
    if (w==NULL)
    {
        fprintf(stderr, "Error -- out of memory\n");
        exit(1);
    }

    sqw = calloc(n, sizeof(double));
    if (sqw==NULL)
    {
        fprintf(stderr, "Error -- out of memory\n");
        exit(1);
    }

    sqwy = calloc(n, sizeof(double));
    if (sqwy==NULL)
    {
        fprintf(stderr, "Error -- out of memory\n");
        exit(1);
    }

    resid = calloc(n, sizeof(double));
    if (resid==NULL)
    {
        fprintf(stderr, "Error -- out of memory\n");
        exit(1);
    }

    /*
     * Initialize quantities before EM iterations
     */

    // Initialize weights for observations
    for (i=0; i<n; i++)
    {
        w[i] = 1;
    }

    /*
     * First iteration
     */

    // Run regression to obtain initial coefficients
    wls(X, n, p, y, w, method, XTX, sqw, sqwX, sqwy, coef);

    // Calculate residuals
    calcResid(X, n, p, y, coef, resid);

    // Calculate tau
    (*tau) = 0;
    for (i=0; i<n; i++)
    {
        (*tau) += resid[i] * resid[i] * w[i];
    }
    (*tau) /= dasum(n, w, 1);

    // Calculate initial log-posterior
    (*logLikelihood) = dt_log(resid, n, nu, 0, sqrt(*tau));

    logLikelihood_tm1 = (*logLikelihood);

    /*
     * EM loop
     */
    for (iter=0; iter<maxIter; iter++)
    {
        // E step: Update u | beta, tau
        for (i=0; i<n; i++)
        {
            w[i] = (nu + 1) / (nu + resid[i]*resid[i]/(*tau));
        }

        // M step: Update beta, tau | u

        // Run regression to obtain coefficients
        wls(X, n, p, y, w, method, XTX, sqw, sqwX, sqwy, coef);

        // Calculate residuals
        calcResid(X, n, p, y, coef, resid);

        // Calculate tau
        (*tau) = 0;
        for (i=0; i<n; i++)
        {
            (*tau) += resid[i] * resid[i] * w[i];
        }
        (*tau) /= dasum(n, w, 1);

        // Calculate log-posterior
        (*logLikelihood) = dt_log(resid, n, nu, 0, sqrt(*tau));

        // Check convergence
        delta = ((*logLikelihood) - logLikelihood_tm1) /
            fabs((*logLikelihood) + logLikelihood_tm1) * 2;

        if (delta < tol)
        {
            logLikelihood_tm1 = (*logLikelihood);
            break;
        }
        logLikelihood_tm1 = (*logLikelihood);
    }

    /*
     * Free allocated memory
     */

    // Free workspace matrices
    free(sqwX);
    sqwX = NULL;

    free(XTX);
    XTX = NULL;

    // Free workspace vectors
    free(w);
    w=NULL;

    free(sqw);
    sqw=NULL;

    free(sqwy);
    sqw=NULL;

    free(resid);
    resid=NULL;

    return iter;
}
