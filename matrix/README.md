matrix
======

Column order matrix implementation suitable for BLAS/LAPACK libraries.

Changes
-------
	New sub-matrix interface that resembles go slices. Creating a sub-matrix does
	not copy elements and changes to submatrix are visible in the original matrix.

	Depreciated row and column handling functions. Use sub-matrix instead.

Todo
----

- ComplexMatrix not up to level of FloatMatrix




