// Copyright (c) Harri Rautila, 2012,2013

// This file is part of github.com/hrautila/matops package. It is free software,
// distributed under the terms of GNU Lesser General Public License Version 3, or
// any later version. See the COPYING tile included in this archive.

package main

import (
	"github.com/henrylee2cn/algorithm/matrix"
	//"github.com/henrylee2cn/algorithm/linalg"
	"flag"
	"fmt"
	"github.com/henrylee2cn/algorithm/linalg/lapack"
	"github.com/henrylee2cn/algorithm/matops"
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

var M, N, KB, MB, NB int
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
var TAUlapack, TAUmatops, Rlapack, Rmatops *matrix.FloatMatrix
var ERRlapack, ERRmatops error

func init() {
	flag.IntVar(&M, "M", 100, "Matrix rows")
	flag.IntVar(&N, "N", 100, "Matrix cols")

	// blocking size; 0 is unblocked versions
	flag.IntVar(&KB, "KB", 0, "Blocking size for blocked invocations")

	// parameters for basic matrix operations
	flag.IntVar(&MB, "MB", 68, "Row blocking size.")
	flag.IntVar(&NB, "NB", 68, "Column blocking size.")
	flag.IntVar(&VPsize, "H", 68, "Viewport size.")
	flag.IntVar(&nWorker, "W", 2, "Number of workers for parallel runs")

	flag.BoolVar(&singleTest, "s", false, "Run single test run for given matrix size.")
	flag.BoolVar(&refTest, "r", false, "Test with lapack reference function.")
	flag.StringVar(&sizeList, "L", "", "Comma separated list of matrix sizes.")
	flag.IntVar(&testCount, "n", 5, "Number of test runs.")

	flag.BoolVar(&check, "C", false, "Check result against lapack reference.")
	flag.BoolVar(&verbose, "v", false, "Be verbose.")
	flag.BoolVar(&asGflops, "g", false, "Report as Gflops.")
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

func checkIPIV(ipiv []int, ipiv2 []int32) bool {
	if len(ipiv) != len(ipiv2) {
		return false
	}
	for k, n := range ipiv {
		if ipiv2[k] != int32(n) {
			return false
		}
	}
	return true
}

// single invocation for matops and lapack functions
func runCheck(A *matrix.FloatMatrix, LB int) (bool, time.Duration, time.Duration) {

	var W *matrix.FloatMatrix = nil
	N := A.Cols()
	tau := matrix.FloatZeros(N, 1)
	if LB > 0 {
		W = matrix.FloatZeros(A.Rows(), LB)
	}
	fnc := func() {
		_, ERRmatops = matops.DecomposeQR(A, tau, W, LB)
	}

	if verbose && N < 10 {
		fmt.Fprintf(os.Stderr, "A start:\n%v\n", A)
	}
	A0 := A.Copy()
	tau0 := tau.Copy()

	mperf.FlushCache()
	time0 := mperf.Timeit(fnc)
	if verbose && N < 10 {
		fmt.Fprintf(os.Stderr, "A end:\n%v\n", A)
		tau.SetSize(1, N, 1)
		fmt.Fprintf(os.Stderr, "tau: %v\n", tau)
	}

	fn2 := func() {
		ERRlapack = lapack.Geqrf(A0, tau0)
	}

	mperf.FlushCache()
	time2 := mperf.Timeit(fn2)
	if verbose && N < 10 {
		fmt.Fprintf(os.Stderr, "A0 end:\n%v\n", A0)
		tau0.SetSize(1, N, 1) // row vector
		fmt.Fprintf(os.Stderr, "tau0: %v\n", tau0)
	}
	// now A == A0 && tau == tau0

	ok := A.AllClose(A0)
	oktau := tau.AllClose(tau0)
	if !ok || !oktau {
		// save result to globals
		Rlapack = A0
		Rmatops = A
		TAUlapack = tau0
		TAUmatops = tau
	}
	return ok && oktau, time0, time2
}

//
func runTest(A *matrix.FloatMatrix, ntest, LB int) time.Duration {
	var W *matrix.FloatMatrix = nil
	var mintime time.Duration

	N := A.Cols()
	tau := matrix.FloatZeros(N, 1)
	if LB > 0 {
		W = matrix.FloatZeros(A.Rows(), LB)
	}
	fnc := func() {
		_, ERRmatops = matops.DecomposeQR(A, tau, W, LB)
	}

	A0 := A.Copy()
	for n := 0; n < ntest; n++ {
		if n > 0 {
			// restore original A
			A0.CopyTo(A)
			tau.Scale(0.0)
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

	var mintime time.Duration

	N := A.Cols()
	tau := matrix.FloatZeros(N, 1)

	fnc := func() {
		ERRlapack = lapack.Geqrf(A, tau)
	}

	A0 := A.Copy()
	for n := 0; n < ntest; n++ {
		if n > 0 {
			// restore original A
			A0.CopyTo(A)
			tau.Scale(0.0)
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
func newMatrix(M, N int) *matrix.FloatMatrix {
	generator := &matrix.NormalFloat{0.0, 2.0}
	A0 := matrix.FloatZeros(M, N)
	A0.SetFrom(generator)
	return A0
}

type testFunc func(*matrix.FloatMatrix, int, int) time.Duration

func runTestsForSizes(test testFunc, sizes []int) map[int]float64 {
	times := make(map[int]float64, len(sizes))
	for _, sz := range sizes {
		A := newMatrix(sz, sz)
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

func columnDiffs(A, B *matrix.FloatMatrix) *matrix.FloatMatrix {
	var c matrix.FloatMatrix
	nrm := matrix.FloatZeros(A.Cols(), 1)
	A0 := A.Copy()
	A0.Minus(B)
	for k := 0; k < A.Cols(); k++ {
		A0.SubMatrix(&c, 0, k, A.Rows(), 1)
		nrm.SetAt(k, 0, matops.Norm2(&c))
	}
	return nrm
}

func rowDiffs(A, B *matrix.FloatMatrix) *matrix.FloatMatrix {
	var r matrix.FloatMatrix
	nrm := matrix.FloatZeros(A.Rows(), 1)
	A0 := A.Copy()
	A0.Minus(B)
	for k := 0; k < A.Rows(); k++ {
		A0.SubMatrix(&r, k, 0, 1, A.Cols())
		nrm.SetAt(k, 0, matops.Norm2(&r))
	}
	return nrm
}

func gFlops(M, N int, secs float64) float64 {
	// flops:
	return 2.0 * float64(int64(M)*int64(N)*int64(N)) / 3.0 / secs * 1e-9
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
		A := newMatrix(M, N)
		Ac := A.Copy()

		if check {
			ok, tm, reftm = runCheck(A, KB)
			if verbose {
				fmt.Fprintf(os.Stderr, "%s: %v\n", testName, tm)
				fmt.Fprintf(os.Stderr, "Reference: [%v] %v (%.2f) \n",
					ok, reftm, tm.Seconds()/reftm.Seconds())
				if !ok {
					fmt.Fprintf(os.Stderr, "ERRlapack: %v\n", ERRlapack)
					fmt.Fprintf(os.Stderr, "ERRmatops: %v\n", ERRmatops)
				}
			}
			if asGflops {
				fmt.Printf("%.4f Gflops [ref: %.4f]\n",
					gFlops(M, N, tm.Seconds()), gFlops(M, N, reftm.Seconds()))
			}
			if !ok {
				//fmt.Printf("A orig:\n%v\n", Ac)
				_ = Ac
				//saveData(Ac)
				cdf := columnDiffs(Rmatops, Rlapack)
				cdf.SetSize(1, cdf.Rows(), 1)
				fmt.Fprintf(os.Stderr, "col diffs: %v\n", cdf)

				rdf := rowDiffs(Rmatops, Rlapack)
				rdf.SetSize(1, rdf.Rows(), 1)
				fmt.Fprintf(os.Stderr, "row diffs: %v\n", rdf)
			}
			return
		}

		if refTest {
			tm = runRefTest(A, testCount, KB)
		} else {
			tm = runTest(A, testCount, KB)
		}

		if asGflops {
			fmt.Printf("%.4f Gflops\n", gFlops(M, N, tm.Seconds()))
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
			times[sz] = gFlops(sz, sz, times[sz])
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
