
int wls(double* X, int n, int p, double* y, double* w, char method,
        double* XTX, double *sqw, double* sqwX, double* sqwy, double* coef);

int calcFitted(double* X, int n, int p,
    double* y,
    double* coef,
    double* fitted);

int calcResid(double* X, int n, int p,
    double* y,
    double* coef,
    double* resid);

