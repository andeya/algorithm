schoolcalc - Support pupils when learning to divide
===================================================

Short description
-----------------
A library for

* dividing integers or rationals. It calculates intermidate steps when performing divisions.
Supports negative numerals.
* Zapfenrechnung. A Zapfenrechnung (Zapfen rechnen) looks like this:

Stand-alone program

* sdivcon: Console program which accepts input from stdin like dividend:divisor and writes the results
  to stdout, displaying the intermediate steps when dividing the pen and paper method
* webzapfen: Web app to display a Zapfenrechnung, including the intermediate division steps 


<pre><code>         27 * 2 = 54
         54 * 3 = 162
        162 * 4 = 648
        648 * 5 = 3240
       3240 * 6 = 19440
      19440 * 7 = 136080
     136080 * 8 = 1088640
    1088640 * 9 = 9797760
    9797760 / 2 = 4898880
    4898880 / 3 = 1632960
    1632960 / 4 = 408240
     408240 / 5 = 81648
      81648 / 6 = 13608
      13608 / 7 = 1944
       1944 / 8 = 243
        243 / 9 = 27</code></pre>

Features
--------

* Division: Arbitrary precision mathematics
* Division: Negative numbers & fractional numbers supported when deviding the pen & paper method
* Division: Intermediate steps are returned for display / control
* ZapfenRechnung: Arbitrary precision mathematics

Installation
------------

Library:

    go get github.com/the42/schoolcalc/schoolcalc

Programs:

Command line program sdivcon to display the intermediate results of SchoolDivide:

      go get github.com/the42/schoolcalc/sdivcon

Web application webzapfen displays a Zapfenrechnung, including the intermediate division steps 

      go get github.com/the42/schoolcalc/webzapfen

Usage
-----

SchoolDivide is covered by test cases. ZapfenRechnung contains a short executable main program
in the comments to the library function.

License
-------

The package is released under the [Simplified BSD
License](http://www.freebsd.org/copyright/freebsd-license.html) See file
"LICENSE"

Testing
-------

To run the tests:

  go test github.com/the42/schoolcalc
