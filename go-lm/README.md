# `go-lm`: Linear models in Go

`go-lm` provides a basic implementation of weighted least squares (WLS) regression and regression with t-distributed residuals. These are implemented using the [`cgo`][cgo] interface, with the C code directly calling standard BLAS/LAPACK functions.

For WLS regression, two methods are provided via the `Wls` function. Setting `method='q'` will use the QR decomposition via the DGELS LAPACK routine. Setting `method='c'` (or anything else) will use the Cholesky decomposition via the DPOSV LAPACK routine.

For linear regression with t-distributed residuals, the optimal PX-EM algorithm of [Meng & van Dyk (1997)][pxem] is implemented via the `LmT` function.

This package requires a [current Go installation][golang] with cgo enabled and the [Atlas BLAS][atlas] and [LAPACK][lapack] libraries to compile. It has been tested on Ubuntu 12.04 with the following packages installed:

```bash
sudo apt-get install libatlas3gf-base libatlas-base-dev liblapack golang build-essential
```

The entire package is provided under the [MIT license][mit].

[cgo]: http://golang.org/cmd/cgo/
[pxem]: http://scholar.google.com/scholar?cluster=8397446410416335201
[atlas]: http://math-atlas.sourceforge.net/
[lapack]: http://www.netlib.org/lapack/
[golang]: http://www.golang.org/
[mit]: http://opensource.org/licenses/MIT
