// Copyright (c) Harri Rautila, 2012,2013

// This file is part of github.com/hrautila/matops package. It is free software,
// distributed under the terms of GNU Lesser General Public License Version 3, or
// any later version. See the COPYING tile included in this archive.

package main

import (
	"flag"
	"fmt"
	"github.com/henrylee2cn/algorithm/matops"
	"github.com/henrylee2cn/algorithm/matrix"
	"github.com/henrylee2cn/algorithm/mperf"
	"io"
	"math/rand"
	"os"
	"runtime"
	"strconv"
	"strings"
	"time"
	//"unsafe"
)

var N, KB, MB, NB int
var randomData bool
var check bool
var verbose bool
var asGflops bool
var asEps bool
var singleTest bool
var testUpper bool
var refTest bool
var noSPD bool
var nWorker int
var testName string
var testCount int
var VPsize int
var sizeList string
var transpose string
var fileName string

// globals for error cases
var Rref, Rmatops *matrix.FloatMatrix
var IPIVref []int
var IPIVmatops []int
var ERRref, ERRmatops error

func init() {
	flag.IntVar(&N, "N", 100, "Matrix rows/cols.")

	// blocking size; 0 is unblocked versions
	flag.IntVar(&KB, "KB", 0, "Blocking size for blocked invocations")

	// parameters for basic matrix operations
	flag.IntVar(&MB, "MB", 68, "Row blocking size.")
	flag.IntVar(&NB, "NB", 68, "Column blocking size.")
	flag.IntVar(&VPsize, "H", 68, "Viewport size.")
	flag.IntVar(&nWorker, "W", 2, "Number of workers for parallel runs")

	flag.BoolVar(&singleTest, "s", false, "Run single test run for given matrix size.")
	flag.BoolVar(&refTest, "r", false, "Test with reference function.")
	flag.StringVar(&sizeList, "L", "", "Comma separated list of matrix sizes.")
	flag.IntVar(&testCount, "n", 5, "Number of test runs.")

	flag.BoolVar(&noSPD, "-nP", false, "Matrix is not SPD.")
	flag.BoolVar(&testUpper, "U", false, "Matrix is UPPER triangular. ")
	flag.BoolVar(&check, "C", false, "Check result against lapack reference.")
	flag.BoolVar(&verbose, "v", false, "Be verbose.")
	flag.BoolVar(&asGflops, "g", false, "Report as Gflops.")
	flag.BoolVar(&asEps, "e", false, "Report as result elements per seconds.")
	flag.StringVar(&testName, "T", "test", "Test name for reporting")
	flag.StringVar(&fileName, "F", "saved.dat", "Filename for source data")
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

func saveData(A *matrix.FloatMatrix) {
	var fd *os.File
	if fileName == "" {
		fileName = testName + ".dat"
	}
	fd, err := os.Create(fileName)
	if err != nil {
		fmt.Fprintf(os.Stderr, "create error: %v\n", err)
		return
	}
	io.WriteString(fd, A.ToString("%14e"))
	fd.Close()
}

func checkIPIV(ipiv []int, ipiv2 []int) bool {
	if len(ipiv) != len(ipiv2) {
		return false
	}
	for k, n := range ipiv {
		if ipiv2[k] != n {
			return false
		}
	}
	return true
}

// single invocation for matops and lapack functions
func runCheck(A *matrix.FloatMatrix, LB int) (bool, time.Duration, time.Duration) {

	var flags matops.Flags
	N := A.Rows()
	ipiv := make([]int, N, N)
	ipiv0 := make([]int, N, N)
	flags = matops.LOWER
	if testUpper {
		flags = matops.UPPER
	}
	W := matrix.FloatZeros(A.Rows(), LB+2)
	fnc := func() {
		_, ERRmatops = matops.DecomposeLDL(A, W, ipiv, flags, LB)
	}

	if verbose && N < 10 {
		fmt.Fprintf(os.Stderr, "A start:\n%v\n", A)
	}
	A0 := A.Copy()
	mperf.FlushCache()
	time0 := mperf.Timeit(fnc)
	if verbose && N < 10 {
		fmt.Fprintf(os.Stderr, "ipiv:%v\n", ipiv)
		fmt.Fprintf(os.Stderr, "A end:\n%v\n", A)
	}

	fn2 := func() {
		_, ERRref = matops.DecomposeLDL(A0, W, ipiv0, flags, 0)
	}
	if verbose && N < 10 {
		fmt.Fprintf(os.Stderr, "A0 start:\n%v\n", A0)
	}
	mperf.FlushCache()
	time2 := mperf.Timeit(fn2)
	if verbose && N < 10 {
		fmt.Fprintf(os.Stderr, "ipiv0:%v\n", ipiv0)
		fmt.Fprintf(os.Stderr, "A0 end:\n%v\n", A0)
	}
	// now A == A0 && ipiv == ipiv0

	ok := A.AllClose(A0)
	okip := checkIPIV(ipiv, ipiv0)
	if !ok || !okip {
		// save result to globals
		Rref = A0
		IPIVref = ipiv0
		Rmatops = A
		IPIVmatops = ipiv
	}
	return ok && okip, time0, time2
}

//
func runTest(A *matrix.FloatMatrix, ntest, LB int) time.Duration {

	var flags matops.Flags
	var mintime time.Duration

	N := A.Rows()
	ipiv := make([]int, N, N)
	flags = matops.LOWER
	if testUpper {
		flags = matops.UPPER
	}

	W := matrix.FloatZeros(A.Rows(), LB+2)
	fnc := func() {
		_, ERRmatops = matops.DecomposeLDL(A, W, ipiv, flags, LB)
	}

	A0 := A.Copy()
	for n := 0; n < ntest; n++ {
		if n > 0 {
			// restore original A
			A0.CopyTo(A)
		}
		mperf.FlushCache()
		time0 := mperf.Timeit(fnc)
		if n == 0 || time0 < mintime {
			mintime = time0
		}
		if verbose {
			fmt.Printf("%.4f ms\n", time0.Seconds()*1000.0)
		}
	}
	return mintime
}

func runRefTest(A *matrix.FloatMatrix, ntest, LB int) time.Duration {

	var flags matops.Flags
	var mintime time.Duration

	N := A.Rows()
	ipiv := make([]int, N, N)
	flags = matops.LOWER
	if testUpper {
		flags = matops.UPPER
	}

	W := matrix.FloatZeros(A.Rows(), LB+2)
	fnc := func() {
		_, ERRref = matops.DecomposeLDL(A, W, ipiv, flags, 0)
	}

	A0 := A.Copy()
	for n := 0; n < ntest; n++ {
		if n > 0 {
			// restore original A
			A0.CopyTo(A)
		}
		mperf.FlushCache()
		time0 := mperf.Timeit(fnc)
		if n == 0 || time0 < mintime {
			mintime = time0
		}
	}
	return mintime
}

// create a new matrix.
func newMatrix(N int) *matrix.FloatMatrix {

	uplo := matrix.Lower
	if testUpper {
		uplo = matrix.Upper
	}
	randNormal := &matrix.NormalFloat{0.0, 2.0}

	A0 := matrix.FloatZeros(N, N)
	if !noSPD {
		A0.SetFrom(randNormal)
		A := matrix.FloatZeros(N, N)
		// A*A.T is positive definite
		matops.Mult(A, A0, A0, 1.0, 0.0, matops.TRANSB)
		return A
	}
	A0.SetFromTrm(generator, uplo)
	return A0
}

type testFunc func(*matrix.FloatMatrix, int, int) time.Duration

func runTestsForSizes(test testFunc, sizes []int) map[int]float64 {

	times := make(map[int]float64, len(sizes))
	for _, sz := range sizes {
		A := newMatrix(sz)
		tm := test(A, testCount, KB)
		times[sz] = tm.Seconds()
	}
	return times
}

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

func gFlops(N int, secs float64) float64 {

	return float64(int64(N)*int64(N)*int64(N)) / (6.0 * secs) * 1e-9
}

func main() {
	flag.Parse()
	runtime.GOMAXPROCS(nWorker)
	matops.NumWorkers(nWorker)
	rand.Seed(time.Now().UnixNano())
	matops.BlockingParams(MB, NB, VPsize)

	var ok bool
	var reftm, tm time.Duration

	if singleTest {
		A := newMatrix(N)
		Ac := A.Copy()

		if check {
			ok, tm, reftm = runCheck(A, KB)
			if verbose {
				fmt.Fprintf(os.Stderr, "%s: %v\n", testName, tm)
				fmt.Fprintf(os.Stderr, "Reference: [%v] %v (%.2f) \n",
					ok, reftm, tm.Seconds()/reftm.Seconds())
				if !ok {
					fmt.Fprintf(os.Stderr, "ERRref:    %v\n", ERRref)
					fmt.Fprintf(os.Stderr, "ERRmatops: %v\n", ERRmatops)
				}
			}
			if asGflops {
				fmt.Printf("%.4f Gflops [ref: %.4f]\n",
					gFlops(N, tm.Seconds()), gFlops(N, reftm.Seconds()))
			}
			if !ok {
				//fmt.Printf("A orig:\n%v\n", Ac)
				saveData(Ac)
			}
			return
		}

		if refTest {
			tm = runRefTest(A, testCount, KB)
		} else {
			tm = runTest(A, testCount, KB)
		}

		if asGflops {
			fmt.Printf("%.4f Gflops\n", gFlops(N, tm.Seconds()))
		} else {
			fmt.Printf("%vs\n", tm.Seconds())
		}
		return
	}

	if len(sizeList) > 0 {
		sizes = parseSizeList(sizeList)
	}
	var times map[int]float64

	if refTest {
		times = runTestsForSizes(runRefTest, sizes)
	} else {
		times = runTestsForSizes(runTest, sizes)
	}
	if asGflops {
		if verbose {
			fmt.Printf("calculating Gflops ...\n")
		}
		for sz := range times {
			times[sz] = gFlops(sz, times[sz])
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
