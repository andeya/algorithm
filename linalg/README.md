This package distributed under LGPL3 license. See file COPYING in repository root.

Linalg package provides thin interface for BLAS/LAPACK libraries. 
It is modelled after CVXOPT python package.  

It is still work in progress and does not implement all LAPACK functions.
And BLAS interfaces for complex valued matrices is not as complete as for float
valued matrices. However, all BLAS functionality is available for complex valued
matrices via generic matrix interface.

It uses real and complex column-major matrix implementation from github.com/hrautila/matrix 

Requirements

external: 
* libblas-dev package   (tested in Ubuntu 12.04)
* liblapack-edv package  (tested in Ubuntu 12.04)

(See the cgo_ files in blas/lapack subdirectories).


other packages:
* go get github.com/hrautila/matrix 




