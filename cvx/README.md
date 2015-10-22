cvx
===

Package cvx is a Go package for solving convex optimization problems.

It is a straightforward translation of parts of the CVXOPT python package for convex
optimization. Spefically it provides interfaces for solving linear and quadratic cone
programs and convex programs with non-linear objectives.

Package cvx depends on column order matrix implementation and access to BLAS and
LAPACK linear algebra libraries. 

To install it do:

* go get github.com/hrautila/matrix
* go get github.com/hrautila/linalg
* go get github.com/hrautila/cvx


For examples see _test.go files. Additional examples and other related material 
see https://github.com/hrautila/go.opt
