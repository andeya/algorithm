PACKAGE

package mahalanobis
    import "github.com/ant0ine/go.mahalanobis"

    Naive implementation of the Mahalanobis distance using go.matrix
    (https://en.wikipedia.org/wiki/Mahalanobis_distance)

    This is me learning Go, it's probably broken, don't use it.

    Example:

	package main
	import (
	    "fmt"
	    "github.com/skelterjohn/go.matrix"
	    "github.com/ant0ine/go.mahalanobis"
	)
	func main() {
	    points := matrix.MakeDenseMatrix([]float64{
	        1, 4, 3, 4,
	        4, 2, 3, 4,
	    }, 2, 4)
	    fmt.Println("4 points:\n", points)
	    target := matrix.MakeDenseMatrix([]float64{
	        3,
	        4,
	    }, 2, 1)
	    fmt.Println("the target point:\n", target)
	    distance, err := mahalanobis.Distance(points, target)
	    if err != nil {
	        panic(err)
	    }
	    fmt.Println("Mahalanobis distance=", distance)
	}

FUNCTIONS

func CovarianceMatrix(points *matrix.DenseMatrix) *matrix.DenseMatrix
    Return the covariance matrix for this set of points (sample covariance
    is used) points.Rows() = dimensions. points.Cols() = number of points.

func Distance(points, target *matrix.DenseMatrix) (float64, error)
    Return the Mahalanobis distance points.Rows() = dimensions.
    points.Cols() = number of points. target.Cols = 1

func DistanceSquare(points, target *matrix.DenseMatrix) (float64, error)
    Return the square of the Mahalanobis distance points.Rows() =
    dimensions. points.Cols() = number of points. target.Cols = 1

func MeanVector(points *matrix.DenseMatrix) *matrix.DenseMatrix
    Given a set a points, return the mean vector. points.Rows() =
    dimensions. points.Cols() = number of points.


