// Naive implementation of the Mahalanobis distance using go.matrix
// (https://en.wikipedia.org/wiki/Mahalanobis_distance)
//
// This is me learning Go, it's probably broken, don't use it.
//
// Example:
//
//      package main
//
//      import (
//          "fmt"
//          "github.com/skelterjohn/go.matrix"
//          "github.com/ant0ine/go.mahalanobis"
//      )
//
//      func main() {
//
//          points := matrix.MakeDenseMatrix([]float64{
//              1, 4, 3, 4,
//              4, 2, 3, 4,
//          }, 2, 4)
//          fmt.Println("4 points:\n", points)
//
//          target := matrix.MakeDenseMatrix([]float64{
//              3,
//              4,
//          }, 2, 1)
//          fmt.Println("the target point:\n", target)
//
//          distance, err := mahalanobis.Distance(points, target)
//          if err != nil {
//              panic(err)
//          }
//          fmt.Println("Mahalanobis distance=", distance)
//      }
package mahalanobis

import (
	"errors"
	"github.com/henrylee2cn/algorithm/go.matrix"
	"math"
)

// Given a set a points, return the mean vector.
// points.Rows() = dimensions. points.Cols() = number of points.
func MeanVector(points *matrix.DenseMatrix) *matrix.DenseMatrix {
	mean := matrix.Zeros(points.Rows(), 1)
	for i := 0; i < points.Rows(); i++ {
		sum := 0.0
		for j := 0; j < points.Cols(); j++ {
			sum += points.Get(i, j)
		}
		mean.Set(i, 0, sum/float64(points.Cols()))
	}
	return mean
}

func sample_covariance_matrix(points, mean *matrix.DenseMatrix) *matrix.DenseMatrix {
	dim := points.Rows()
	cov := matrix.Zeros(dim, dim)
	for i := 0; i < dim; i++ {
		for j := 0; j < dim; j++ {
			if i > j {
				// symetric matrix
				continue
			}
			// TODO in go routines ?
			sum := 0.0
			for k := 0; k < points.Cols(); k++ {
				sum += (points.Get(i, k) - mean.Get(i, 0)) * (points.Get(j, k) - mean.Get(j, 0))
			}

			// this is the sample covariance, divide by (N - 1)
			covariance := sum / (float64(points.Cols() - 1))

			cov.Set(i, j, covariance)
			// symetric matrix
			cov.Set(j, i, covariance)

		}
	}
	return cov
}

// Return the covariance matrix for this set of points (sample covariance is used)
// points.Rows() = dimensions. points.Cols() = number of points.
func CovarianceMatrix(points *matrix.DenseMatrix) *matrix.DenseMatrix {
	mean := MeanVector(points)
	return sample_covariance_matrix(points, mean)
}

// Return the square of the Mahalanobis distance
// points.Rows() = dimensions. points.Cols() = number of points.
// target.Cols = 1
func DistanceSquare(points, target *matrix.DenseMatrix) (float64, error) {
	// TODO support multiple points for target, and return a matrix of distances
	if target.Rows() != points.Rows() {
		err := errors.New("target does not have the same dimension than points")
		return -1, err
	}

	mean := MeanVector(points)

	delta := target.Copy()
	delta.SubtractDense(mean)

	cov := sample_covariance_matrix(points, mean)

	inv, err := cov.Inverse()
	if err != nil {
		return -1, err
	}

	product1, err := inv.TimesDense(delta)
	if err != nil {
		return -1, err
	}
	delta_t := delta.Transpose()
	product2, err := delta_t.TimesDense(product1)
	if err != nil {
		return -1, err
	}

	return product2.Get(0, 0), nil
}

// Return the Mahalanobis distance
// points.Rows() = dimensions. points.Cols() = number of points.
// target.Cols = 1
func Distance(points, target *matrix.DenseMatrix) (float64, error) {
	square, err := DistanceSquare(points, target)
	if err != nil {
		return -1, err
	}
	return math.Sqrt(square), nil

}
