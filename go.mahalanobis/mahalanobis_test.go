package mahalanobis

import (
	//    "fmt"
	"github.com/henrylee2cn/algorithm/go.matrix"
	"math"
	"testing"
)

func is_near(value, expected, epsilon float64) bool {
	return math.Abs(value-expected) < epsilon
}

func TestMeanVector(t *testing.T) {

	points := matrix.MakeDenseMatrix([]float64{
		1, 1, 1,
		1, 1, 1,
	}, 2, 3)
	expected := matrix.MakeDenseMatrix([]float64{
		1,
		1,
	}, 2, 1)
	result := MeanVector(points)
	if !matrix.Equals(result, expected) {
		t.Error()
	}

	points = matrix.MakeDenseMatrix([]float64{
		0, 1, 2,
		0, 2, 4,
	}, 2, 3)
	expected = matrix.MakeDenseMatrix([]float64{
		1,
		2,
	}, 2, 1)
	result = MeanVector(points)
	if !matrix.Equals(result, expected) {
		t.Error()
	}
}

func TestCovarianceMatrix(t *testing.T) {

	// no (co)variance
	// R: var(cbind(c(1, 1), c(1, 1)))
	points := matrix.MakeDenseMatrix([]float64{
		1, 1,
		1, 1,
	}, 2, 2)
	expected := matrix.MakeDenseMatrix([]float64{
		0, 0,
		0, 0,
	}, 2, 2)
	result := CovarianceMatrix(points)
	//fmt.Println("covariance:\n", result)
	if !matrix.Equals(result, expected) {
		t.Error()
	}

	// diagonale case
	// R: var(cbind(c(0, 4, 2, 2), c(2, 2, 0, 4)))
	points = matrix.MakeDenseMatrix([]float64{
		0, 4, 2, 2,
		2, 2, 0, 4,
	}, 2, 4)
	expected = matrix.MakeDenseMatrix([]float64{
		2.66, 0,
		0, 2.66,
	}, 2, 2)
	result = CovarianceMatrix(points)
	//fmt.Println("covariance:\n", result)
	if !matrix.ApproxEquals(result, expected, 0.01) {
		t.Error()
	}

	// another case
	// R: var(cbind(c(9, 3, 5), c(3, 4, 1)))
	points = matrix.MakeDenseMatrix([]float64{
		9, 3, 5,
		3, 4, 1,
	}, 2, 3)
	expected = matrix.MakeDenseMatrix([]float64{
		9.33, -0.66,
		-0.66, 2.33,
	}, 2, 2)
	result = CovarianceMatrix(points)
	//fmt.Println("covariance:\n", result)
	if !matrix.ApproxEquals(result, expected, 0.01) {
		t.Error()
	}
}

func TestDistance(t *testing.T) {

	// R: x = cbind(c(9, 3, 5), c(3, 4, 1))
	points := matrix.MakeDenseMatrix([]float64{
		9, 3, 5,
		3, 4, 1,
	}, 2, 3)
	target := matrix.MakeDenseMatrix([]float64{
		1,
		1,
	}, 2, 1)

	square, err := DistanceSquare(points, target)
	if err != nil {
		t.Fatal(err)
	}
	var expected float64
	expected = 4.08
	if !is_near(square, expected, 0.01) {
		t.Error()
	}

	// R: mahalanobis(c(1,1), colMeans(x), var(x))
	distance, err := Distance(points, target)
	if err != nil {
		t.Fatal(err)
	}
	expected = 2.02
	if !is_near(distance, expected, 0.01) {
		t.Error()
	}
}
