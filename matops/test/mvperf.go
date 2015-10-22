// Copyright (c) Harri Rautila, 2012,2013

// This file is part of github.com/hrautila/matops package. It is free software,
// distributed under the terms of GNU Lesser General Public License Version 3, or
// any later version. See the COPYING tile included in this archive.

package main

import (
	"github.com/henrylee2cn/algorithm/linalg"
	"github.com/henrylee2cn/algorithm/linalg/blas"
	"github.com/henrylee2cn/algorithm/matrix"
	"github.com/henrylee2cn/algorithm/mperf"
	//"github.com/henrylee2cn/algorithm/matops/calgo"
	"flag"
	"fmt"
	"github.com/henrylee2cn/algorithm/matops"
	"math/rand"
	"os"
	"runtime"
	"strconv"
	"strings"
	"time"
	//"unsafe"
)

var M, N, P, MB, NB int
var randomData bool
var check bool
var verbose bool
var asGflops bool
var asEps bool
var singleTest bool
var nWorker int
var testName string
var testCount int
var VPsize int
var sizeList string
var transpose string

func init() {
	flag.IntVar(&M, "M", 600, "Matrix A rows.")
	flag.IntVar(&N, "N", 600, "Matrix B cols.")
	flag.IntVar(&P, "P", 600, "Matrix A cols, B rows.")
	flag.IntVar(&MB, "MB", 64, "Row blocking size.")
	flag.IntVar(&NB, "NB", 64, "Column blocking size.")
	flag.IntVar(&nWorker, "W", 2, "Number of workers for parallel runs")
	flag.IntVar(&VPsize, "H", 64, "Viewport size.")
	flag.BoolVar(&check, "C", false, "Check result against reference (gemm).")
	flag.BoolVar(&verbose, "v", false, "Be verbose.")
	flag.BoolVar(&asGflops, "g", false, "Report as Gflops.")
	flag.BoolVar(&asEps, "e", false, "Report as result elements per second.")
	flag.BoolVar(&randomData, "R", true, "Generate random data.")
	flag.BoolVar(&singleTest, "s", false, "Run single test run for given matrix sizes.")
	flag.IntVar(&testCount, "n", 5, "Number of test runs.")
	flag.StringVar(&testName, "T", "test", "Test name for reporting")
	flag.StringVar(&sizeList, "L", "", "Comma separated list of sizes.")
	flag.StringVar(&transpose, "t", "N", "Transpose op, N, A, B, AB")
}

var sizes []int = []int{
	10, 30, 50, 70, 90,
	100, 200, 300, 400, 500, 600, 700, 800, 900,
	1000, 1100, 1200, 1300, 1400, 1500}

func index(i, r, sz int) int {
	if i == r {
		return sz
	}
	return i*sz/r - ((i * sz / r) & 0x3)
}

func TestTemplate(m, n, p int) (fnc func(), A, X, Y *matrix.FloatMatrix) {
	A = matrix.FloatNormal(m, n)
	X = matrix.FloatNormal(n, 1)
	Y = matrix.FloatZeros(m, 1)
	fnc = func() {
		// test core here
	}
	return
}

func CTestMVMult(m, n, p int) (fnc func(), A, X, Y *matrix.FloatMatrix) {
	A = matrix.FloatNormal(m, n)
	X = matrix.FloatNormal(n, 1)
	Y = matrix.FloatZeros(m, 1)
	fnc = func() {
		matops.MVMult(Y, A, X, 1.0, 1.0)
	}
	return
}

func CTestMVMultTransA(m, n, p int) (fnc func(), A, X, Y *matrix.FloatMatrix) {
	A = matrix.FloatNormal(n, m)
	X = matrix.FloatNormal(n, 1)
	Y = matrix.FloatZeros(m, 1)
	fnc = func() {
		matops.MVMultTransA(Y, A, X, 1.0, 1.0)
	}
	return
}

func CTestGemv(m, n, p int) (fnc func(), A, X, Y *matrix.FloatMatrix) {
	A = matrix.FloatNormal(m, n)
	X = matrix.FloatNormal(n, 1)
	Y = matrix.FloatZeros(m, 1)
	fnc = func() {
		blas.GemvFloat(A, X, Y, 1.0, 1.0)
	}
	return
}

func CTestGemvTransA(m, n, p int) (fnc func(), A, X, Y *matrix.FloatMatrix) {
	A = matrix.FloatNormal(n, m)
	X = matrix.FloatNormal(n, 1)
	Y = matrix.FloatZeros(m, 1)
	A = A.Transpose()
	fnc = func() {
		blas.GemvFloat(A, X, Y, 1.0, 1.0, linalg.OptTrans)
	}
	return
}

func CheckNoTrans(A, X, Y *matrix.FloatMatrix) {
	blas.GemvFloat(A, X, Y, 1.0, 1.0)
}

func CheckTransA(A, X, Y *matrix.FloatMatrix) {
	blas.GemvFloat(A, X, Y, 1.0, 1.0, linalg.OptTrans)
}

var tests map[string]mperf.MatrixTestFunc = map[string]mperf.MatrixTestFunc{
	"MVMult":       CTestMVMult,
	"MVMultTransA": CTestMVMultTransA,
	"GemvTransA":   CTestGemvTransA,
	"Gemv":         CTestGemv}

func parseSizeList(s string) []int {
	sl := strings.Split(s, ",")
	il := make([]int, 0)
	for _, snum := range sl {
		n, err := strconv.ParseInt(snum, 0, 32)
		if err == nil {
			il = append(il, int(n))
		}
	}
	return il
}

func main() {
	flag.Parse()
	runtime.GOMAXPROCS(nWorker)
	matops.NumWorkers(nWorker)
	rand.Seed(time.Now().UnixNano())

	testFunc, ok := tests[testName]
	if !ok {
		fmt.Printf("Error: test %s does not exists.\nKnown tests:\n", testName)
		for tname := range tests {
			fmt.Printf("\t%s\n", tname)
		}
		return
	}
	var checkFunc mperf.MatrixCheckFunc
	if transpose[0] == 'A' {
		checkFunc = CheckTransA
	} else {
		checkFunc = CheckNoTrans
	}

	if singleTest {
		fnc, A, X, Y0 := testFunc(M, N, P)
		mperf.FlushCache()
		tm := mperf.Timeit(fnc)
		if check {
			reftime, ok := mperf.CheckWithFunc(A, X, Y0, checkFunc)
			if verbose {
				fmt.Fprintf(os.Stderr, "%s: %v\n", testName, tm)
				fmt.Fprintf(os.Stderr, "Reference: [%v] %v (%.2f) \n",
					ok, reftime, tm.Seconds()/reftime.Seconds())
			}
		}
		//sec, _ := mperf.SingleTest(testName, testFunc, M, N, P, check, verbose)
		if asGflops {
			gflops := 2.0 * float64(int64(M)*int64(N)*int64(P)) / tm.Seconds() * 1e-9
			fmt.Printf("%.4f Gflops\n", gflops)
		} else if asEps {
			eps := float64(int64(M)*int64(N)) / tm.Seconds()
			fmt.Printf("%.4f Eps\n", eps)
		} else {
			fmt.Printf("%vs\n", tm.Seconds())
		}
		return
	}

	if len(sizeList) > 0 {
		sizes = parseSizeList(sizeList)
	}
	times := mperf.MultipleSizeTests(testFunc, sizes, testCount, verbose)
	if asGflops || asEps {
		if verbose {
			fmt.Printf("calculating Gflops ...\n")
		}
		for sz := range times {
			n := int64(sz)
			if asGflops {
				times[sz] = 2.0 * float64(n*n) / times[sz] * 1e-9
			} else {
				times[sz] = float64(n) / times[sz]
			}
		}
	}
	// print out as python dictionary
	fmt.Printf("{")
	i := 0
	for sz := range times {
		if i > 0 {
			fmt.Printf(", ")
		}
		fmt.Printf("%d: %v", sz, times[sz])
		i++
	}
	fmt.Printf("}\n")
}

// Local Variables:
// tab-width: 4
// indent-tabs-mode: nil
// End:
