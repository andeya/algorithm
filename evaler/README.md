evaler
======

[![Build Status](https://travis-ci.org/soniah/evaler.svg?branch=master)](https://travis-ci.org/soniah/evaler)
[![Coverage](http://gocover.io/_badge/github.com/soniah/evaler)](http://gocover.io/github.com/soniah/evaler)
[![GoDoc](https://godoc.org/github.com/soniah/evaler?status.png)](http://godoc.org/github.com/soniah/evaler)
https://github.com/soniah/evaler

Package evaler implements a simple floating point arithmetic expression evaluator.

Evaler uses Dijkstra's Shunting Yard algorithm [1] to convert an infix
expression to postfix/RPN format [2], then evaluates the RPN expression. The
implementation is adapted from a Java implementation at [3]. The results are
returned as a `math/big *big.Rat`.

This is release 2.0. The previous version that returned results as float64
is in the branch float64.

Usage
-----

```go
result, err := evaler.Eval("1+2")
```

Operators
---------

The operators supported are:

```+ - * / ** () < >```

< (less than) and > (greater than) will get lowest precedence, all
other precedence is as expected (BODMAS [4]).

< and > tests will evaluate to 0.0 for false and 1.0 for true, allowing
expressions like:

```
3 * (1 < 2) # returns 3.0
3 * (1 > 2) # returns 0.0
```

Minus implements both binary and unary operations (thanks @hiroxy).

Issues
------

The math/big library doesn't have an exponent function (**), and implenting one
for big.Rat numbers is non-trivial. As a work around, arguments are converted
to float64's, the calculation is done using the math.Pow() function, the
result is converted to a big.Rat and placed back on the stack.

Documentation
-------------

http://godoc.org/github.com/soniah/evaler

There are also a number of utility functions (eg BigratToFloat(),
BigratToInt()) that may be useful when working with evaler.

Author
------

Sonia Hamilton

http://blog.snowfrog.net

sonia@snowfrog.net

License
-------

Modified BSD License (BSD-3)

Links
-----

[1] http://en.wikipedia.org/wiki/Shunting-yard_algorithm

[2] http://en.wikipedia.org/wiki/Reverse_Polish_notation

[3] http://willcode4beer.com/design.jsp?set=evalInfix

[4] http://www.mathsisfun.com/operation-order-bodmas.html

