// Copyright (c) Harri Rautila, 2012

// This file is part of github.com/hrautila/cvx package.
// It is free software, distributed under the terms of GNU Lesser General Public
// License Version 3, or any later version. See the COPYING tile included in this archive.

// Package for saving and check execution variables at checkpoints.
package checkpnt

import (
	"bufio"
	"errors"
	"fmt"
	"github.com/henrylee2cn/algorithm/cvx/sets"
	"github.com/henrylee2cn/algorithm/linalg/blas"
	"github.com/henrylee2cn/algorithm/matrix"
	"math"
	"os"
	"strconv"
	"strings"
)

type Verifiability interface {
	// Verify variable against reference values. Returns difference.
	Verify(vals ...interface{}) float64
}
type Verifiable interface {
	// Verify variable against reference values. Returns difference.
	Verify(dataline string) float64
	ShowError(dataline string)
}

type checkpoint struct {
	name     string
	filepath string
	major    int
	minor    int
}

type dataPoint struct {
	mtx      *matrix.FloatMatrix
	w        *sets.FloatMatrixSet
	fvar     *float64
	vvar     Verifiable
	panicVar bool
	inErrorM bool
	inErrorF bool
	ckp      *checkpoint
}

type variableTable map[string]*dataPoint
type dataTable map[string]string

var variables variableTable
var active bool = false
var spmajor int
var sppath string
var normError float64
var diffError float64
var verbose bool
var minorstack []int
var minorpointer int
var spformat string

func init() {
	variables = make(variableTable, 20)
	spmajor = 0
	sppath = "./"
	normError = 1e-15
	diffError = 1e-12
	verbose = false
	active = false
	minorstack = make([]int, 20)
	minorpointer = 0
	spformat = "%12.5e"
}

type mVariable struct {
	mtx *matrix.FloatMatrix
}

func (m *mVariable) Verify(dataline string) float64 {
	refdata, err := matrix.FloatParse(dataline)
	if err != nil {
		fmt.Printf("parse error: %s", err)
		return 0.0
	}
	return blas.Nrm2Float(matrix.Minus(m.mtx, refdata))
}

func (m *mVariable) ShowError(dataline string) {
	refdata, err := matrix.FloatParse(dataline)
	if err != nil {
		fmt.Printf("parse error: %s", err)
		return
	}
	df := matrix.Minus(m.mtx, refdata)
	emtx, _ := matrix.FloatMatrixStacked(matrix.StackRight, m.mtx, refdata, df)
	fmt.Printf("my data | ref.data | diff \n%v\n", emtx.ToString(spformat))
}

type fVariable struct {
	fvar *float64
}

func (f *fVariable) Verify(dataline string) float64 {
	refdata, err := strconv.ParseFloat(strings.Trim(dataline, " "), 64)
	if err != nil {
		fmt.Printf("parse error: %s", err)
		return 0.0
	}
	return math.Sqrt((*f.fvar - refdata) * (*f.fvar - refdata))
}

func (f *fVariable) ShowError(dataline string) {
	refdata, err := strconv.ParseFloat(strings.Trim(dataline, " "), 64)
	if err != nil {
		fmt.Printf("parse error: %s", err)
		return
	}
	df := (*f.fvar - refdata)
	s := fmt.Sprintf("my data | ref.data | diff \n %s %s %s\n", spformat, spformat, spformat)
	fmt.Printf(s, *f.fvar, refdata, df)
}

// Return current major number.
func Major() int {
	return spmajor
}

// Advance major number by one.
func MajorNext() {
	if active {
		spmajor += 1
	}
}

// Push new minor number on to stack.
func MinorPush(minor int) {
	if !active {
		return
	}
	if minorpointer == len(minorstack) {
		// stack full
		return
	}
	minorstack[minorpointer] = minor
	minorpointer += 1
}

// Pop minor number on top of the stack.
func MinorPop() int {
	if !active {
		return -1
	}
	if minorpointer == 0 {
		// stack empty
		return -1
	}
	minorpointer -= 1
	return minorstack[minorpointer]
}

// Get minor number on top of the stack.
func MinorTop() int {
	if !active {
		return -1
	}
	if minorpointer == 0 {
		// stack empty,
		return -1
	}
	return minorstack[minorpointer-1]
}

// Test if minor number stack is empty.
func MinorEmpty() bool {
	return minorpointer == 0
}

// Add matrix variable as checkpointable variable.
func AddMatrixVar(name string, mtx *matrix.FloatMatrix) {
	if !active {
		return
	}
	_, ok := variables[name]
	if !ok {
		variables[name] = &dataPoint{vvar: &mVariable{mtx}}
	}
}

// Set or unset panic flag for variable.
func PanicVar(name string, ispanic bool) {
	if !active {
		return
	}
	_, ok := variables[name]
	if ok {
		variables[name].panicVar = ispanic
	}
}

// Add or update float variable as check point variable.
func AddFloatVar(name string, fptr *float64) {
	if !active {
		return
	}
	variables[name] = &dataPoint{vvar: &fVariable{fptr}}
}

// Add or update float variable as check point variable.
func AddVerifiable(name string, vvar Verifiable) {
	if !active {
		return
	}
	variables[name] = &dataPoint{vvar: vvar}
}

// Add or update scaling matrix set to checkpoint variables.
func AddScaleVar(w *sets.FloatMatrixSet) {
	if !active {
		return
	}
	// add all matrices of scale set to variable table
	for _, key := range w.Keys() {
		mset := w.At(key)
		for k, m := range mset {
			name := fmt.Sprintf("%s.%d", key, k)
			variables[name] = &dataPoint{vvar: &mVariable{m}}
		}
	}
}

// Print checkpoint variables.
func PrintVariables() {
	for name := range variables {
		dp := variables[name]
		if dp.mtx != nil {
			fmt.Printf("'%s' matrix (%d, %d)\n", name, dp.mtx.Rows(), dp.mtx.Cols())
		} else if dp.w != nil {
			fmt.Printf("'%s' matrix set \n", name)
		}
	}
}

// Report on check point variables. Prints out the last check point variable turned invalid.
func Report() {
	if !active {
		return
	}
	for name := range variables {
		dp := variables[name]
		if dp.inErrorM {
			fmt.Printf("%8s invalidated at %d.%04d %-12s [%s]\n",
				name, dp.ckp.major, dp.ckp.minor, dp.ckp.name, dp.ckp.filepath)
		}
		if dp.inErrorF {
			fname := name
			if dp.mtx != nil && dp.fvar != nil {
				fname = name + ".t"
			}
			fmt.Printf("%8s invalidated at %d.%04d %-12s [%s]\n",
				fname, dp.ckp.major, dp.ckp.minor, dp.ckp.name, dp.ckp.filepath)
		}
	}
}

func Format(format string) {
	spformat = format
}

func Reset(path string) {
	for name := range variables {
		delete(variables, name)
	}
	spmajor = 0
	active = false
	sppath = path
}

func Activate() {
	active = true
}

func Verbose(flag bool) {
	verbose = flag
}

// Check variables at checkpoint.
func Check(name string, minor int) {
	if !active {
		return
	}
	data, checkp, err := readCkpData(name, minor)
	if err != nil {
		//fmt.Printf("error when reading savepoint '%s': %v\n", name, err)
		return
	}
	// loop through all our savepoint variables
	for varname := range variables {
		dataline, ok := (*data)[varname]
		if !ok {
			// varname not found in savepoint reference data
			continue
		}
		// internal data value
		mydp := variables[varname]
		if mydp.vvar == nil {
			continue
		}
		diff := mydp.vvar.Verify(dataline)
		if diff > diffError {
			if !mydp.inErrorM {
				fmt.Printf("%d.%d sp '%s'[file:%s] variable '%s': normError = %9.2e\n",
					checkp.major, checkp.minor, checkp.name, checkp.filepath, varname, diff)
				mydp.ckp = checkp
			}
			if verbose && !mydp.inErrorM || mydp.panicVar {
				mydp.vvar.ShowError(dataline)
			}
			mydp.inErrorM = true
			if mydp.panicVar {
				panic("variable divergence error ...")
			}
		} else {
			if mydp.inErrorM {
				fmt.Printf("%d.%d sp '%s'[file:%s] variable '%s': returned to valid\n",
					checkp.major, checkp.minor, checkp.name, checkp.filepath, varname)
				mydp.ckp = nil
			}
			mydp.inErrorM = false
		}
	}
}

func readCkpData(name string, minor int) (data *dataTable, ckp *checkpoint, err error) {
	err = nil
	refname := ""
	path := fmt.Sprintf("%s/%04d-%04d.%s", sppath, spmajor, minor, name)
	file, ferr := os.Open(path)
	if ferr != nil {
		data = nil
		err = ferr
		return
	}

	refdata := make(dataTable, 30)
	linereader := bufio.NewReader(file)
	reading := true
	lineno := 0
	var lerr error = nil
	for reading {
		line := ""
		isPrefix := true
		for isPrefix {
			var data []byte
			data, isPrefix, lerr = linereader.ReadLine()
			if lerr != nil {
				reading = false
				goto ready
			}
			line += string(data)
		}
		if lineno == 0 {
			index := strings.Index(line, ":")
			refname = strings.Trim(line[index+1:], " \n")
			if refname != name {
				err = errors.New(fmt.Sprintf("expecting sp: %s, found %s", name, refname))
				return
			}
		} else {
			vartype := parseLine(line, &refdata)
			_ = vartype
		}
		lineno += 1
	}
ready:
	data = &refdata
	ckp = &checkpoint{name: name, filepath: path, major: spmajor, minor: minor}
	return
}

func parseLine(line string, data *dataTable) string {
	index := strings.Index(line, ":")
	if index == -1 {
		// unknown line
		return ""
	}
	prefix := strings.Trim(line[:index], " ")
	pparts := strings.Fields(prefix)
	if len(pparts) != 3 {
		// corrupted prefix
		return ""
	}
	(*data)[pparts[0]] = strings.Trim(line[index+1:], " ")
	return pparts[1]
}

// Local Variables:
// tab-width: 4
// End:
