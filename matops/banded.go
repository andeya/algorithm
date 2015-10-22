package matops

/*

UPPER TRIANGULAR-BAND-PACKED STORAGE

Following is an example of an upper triangular band matrix U of order 6 and an upper
band width of 3. (M=6, N=6, K=3)

Given the following matrix U:

                    | 11  12  13  14   0   0 |
                    |  0  22  23  24  25   0 |
                    |  0   0  33  34  35  36 |
                    |  0   0   0  44  45  46 |
                    |  0   0   0   0  55  56 |
                    |  0   0   0   0   0  66 |

you store it in upper-triangular-band-packed storage mode in array UTB, declared as UTB(4,6),
as follows:


                    |  *   *   *  14  25  36 |
             UTB =  |  *   *  13  24  35  46 |
                    |  *  12  23  34  45  56 |
                    | 11  22  33  44  55  66 |

LOWER TRIANGULAR-BAND-PACKED STORAGE

Following is an example of a lower triangular band matrix L of order 6 and a lower
band width of 2. (M=6, N=6, K=2)

Given the following matrix L:

                    | 11   0   0   0   0   0 |
                    | 21  22   0   0   0   0 |
                    | 31  32  33   0   0   0 |
                    |  0  42  43  44   0   0 |
                    |  0   0  53  54  55   0 |
                    |  0   0   0  64  65  66 |

you store it in lower-triangular-band-packed storage mode in array LTB, declared as LTB(3,6),
as follows:

                    | 11  22  33  44  55  66 |
             LTB =  | 21  32  43  54  65   * |
                    | 31  42  53  64   *   * |


*/

import (
	"github.com/henrylee2cn/algorithm/matrix"
)

// Convert triangular band matrix S to new matrix R with elements stored in
// triangular-band-packed mode. Returned matrix has dimensions R.Rows() == K+1
// and R.Cols() == S.Cols(). Parameter flags must have either UPPER or LOWER bit set.
func BandedTrmMatrix(S *matrix.FloatMatrix, K int, flags Flags) (R *matrix.FloatMatrix) {
	if S.Rows() != S.Cols() {
		return nil
	}
	M := S.Rows()
	N := S.Cols()

	R = nil
	Sr := S.FloatArray()
	if flags&UPPER != 0 {
		// Upper triangular matrix
		R = matrix.FloatZeros(K+1, M)
		Rr := R.FloatArray()
		for j := 0; j < N; j++ {
			m := K + 1 - j
			// is = max(0, j-K)
			is := j - K
			if is < 0 {
				is = 0
			}
			for i := is; i <= j; i++ {
				Rr[j*(K+1)+m+i-1] = Sr[j*M+i]
			}
		}
	} else if flags&LOWER != 0 {
		// Lower triangular matrix
		R = matrix.FloatZeros(K+1, M)
		Rr := R.FloatArray()
		for j := 0; j < N; j++ {
			m := 1 - j
			// ie = min(N, j+K+1)
			ie := j + K + 1
			if ie >= N {
				ie = N
			}
			for i := j; i < ie; i++ {
				Rr[j*(K+1)+m+i-1] = Sr[j*M+i]
			}
		}
	}
	return
}

// Local Variables:
// tab-width: 4
// indent-tabs-mode: nil
// End:
