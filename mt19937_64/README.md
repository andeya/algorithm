mt19937_64
==========

Mersenne Twister (int64) for Go

Usage
-----
```go
import "github.com/farces/mt19937_64"
...
x := mt19937_64.New()
(optional) x.Seed(int64)
result := x.Generate()
```

or, if being used as a rand.Source:

```go
x := rand.New(mt19937_64.New())
... any rand methods are now available using mt_64 ...
```

Initial seed value (if no seed provided) is 5489.

Based on the C Implementation by Makoto Matsumoto and Takuji Nishimura (c) 2004
see: http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/VERSIONS/C-LANG/mt19937-64.c