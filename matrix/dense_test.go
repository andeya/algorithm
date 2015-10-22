package matrix

import (
    //"fmt"
    "testing"
)

var X, A, B *FloatMatrix

func TestFParse(t *testing.T) {
    t.Logf("Test matrix string parsing (MatLab style).\n")
    s := `[1.0 2.0 3.0; 4.0 5.0 6.0]`
    A, err := FloatParse(s)
    if err != nil {
        t.Fatal(err)
    }
    t.Logf("A :\n%v\n", A)
}

func TestFParse2(t *testing.T) {
    t.Logf("Test matrix string parsing (Python style).\n")
    s2 := "[-7.44e-01  1.11e-01  1.29e+00  2.62e+00 -1.82e+00]" +
        "[ 4.59e-01  7.06e-01  3.16e-01 -1.06e-01  7.80e-01]" +
        "[-2.95e-02 -2.22e-01 -2.07e-01 -9.11e-01 -3.92e-01]" +
        "[-7.75e-01  1.03e-01 -1.22e+00 -5.74e-01 -3.32e-01]" +
        "[-1.80e+00  1.24e+00 -2.61e+00 -9.31e-01 -6.38e-01]"

    A, err := FloatParse(s2)
    if err != nil {
        t.Fatal(err)
    }
    t.Logf("Py-A :\n%v\n", A)
    // this produces error (column count mismatch)
    s := "[1.0  2.0  3.0 4.0]\n[1.1  2.1 3.1]"
    A, err = FloatParse(s)
    if err != nil {
        t.Logf("error: %v\n", err)
    }
}

func TestCRandom(t *testing.T) {
    B := FloatUniform(3, 2)
    A := FloatUniformSymmetric(3)
    t.Logf("B:\n%v\n", B)
    t.Logf("A symmetric:\n%v\n", A)
}

func TestCCopy(t *testing.T) {
    t.Logf("Test creating and setting elements.\n")
    A := FloatNew(2, 3, []float64{1, 4, 2, 5, 3, 6})
    t.Logf("A:\n%v\n", A)
    C := A.Copy()
    t.Logf("C:\n%v\n", C)
    C.SetAt(0, 1, 10.0*C.GetAt(0, 1))
    B := FloatNew(3, 2, []float64{1, 2, 3, 4, 5, 6})
    t.Logf("B:\n%v\n", B)
}

func TestFBool(t *testing.T) {
    t.Logf("Test matrix boolean operators.\n")
    A := FloatNew(3, 2, []float64{1, 4, 2, 5, 3, 6})
    t.Logf("A:\n%v\n", A)
    B := A.Copy()
    //B := MakeMatrix(3, 2, []float64{1,2,3,4,5,6})
    t.Logf("B:\n%v\n", B)
    t.Logf("A == B: %v\n", A.Equal(B))
    t.Logf("A <  B: %v\n", A.Less(B))
    t.Logf("A <= B: %v\n", A.LessOrEqual(B))
    t.Logf("A >  B: %v\n", A.Greater(B))
    t.Logf("A >= B: %v\n", A.GreaterOrEqual(B))
}

func TestFMath(t *testing.T) {
    t.Logf("Test matrix basic math.\n")
    A := FloatZeros(2, 2)
    t.Logf("A\n%v\n", A)
    A.Add(1.0)
    t.Logf("A += 1.0\n%v\n", A)
    A.Scale(9.0)
    t.Logf("A *= 9.0\n%v\n", A)
    A.Add(-1.0)
    t.Logf("A -= 1.0\n%v\n", A)
}

func TestMath2(t *testing.T) {
    m := FloatOnes(8, 1)
    m.Add(1.0, 0, 1, 2)
    m.Add(5.0, 5, 6, 7)
    t.Logf("%v\n", m)
}

func TestTimes(t *testing.T) {
    A := FloatWithValue(4, 4, 1.0)
    x := FloatWithValue(4, 1, 2.0)
    t.Logf("A*x=\n%v\n", A.Times(x))
    B := FloatWithValue(4, 4, 1.0)
    t.Logf("A*B=\n%v\n", A.Times(B))
}

func TestBigVector(t *testing.T) {
    A := FloatDiagonal(20000, 1.0)
    X := FloatWithValue(20000, 1, 2.0)
    t.Logf("A*x=%d elements\n", Times(A, X).NumElements())
}

func TestTimesBig(t *testing.T) {
    A := FloatDiagonal(800, 1.0)
    B := FloatWithValue(800, 800, 2.0)
    C := Times(A, B)
    t.Logf("A*B = %d elements\n", C.NumElements())
}

func TestTimesBig2(t *testing.T) {
    A := FloatWithValue(6000, 10, 1.0)
    B := FloatWithValue(10, 6000, 1.0)
    C := Times(A, B)
    t.Logf("A*B = %d elements\n", C.NumElements())
}

func TestTimesBig3(t *testing.T) {
    A := FloatWithValue(6000, 10, 1.0)
    B := FloatWithValue(10, 6000, 1.0)
    C := Times(B, A)
    t.Logf("B*A = %d elements\n", C.NumElements())
}

func TestFuncs(t *testing.T) {
    t.Logf("Test matrix element wise math.\n")
    A := FloatZeros(2, 3)
    AddTwo := func(n float64) float64 {
        return n + 2.0
    }
    C := Apply(A, AddTwo)
    t.Logf("C = AddTwo(A):\n%v\n", C)
}

func TestIndexing(t *testing.T) {
    A := FloatVector([]float64{0, 1, 2, 3, 4, 5})
    t.Logf(" 0: %v\n", A.GetIndex(0))
    t.Logf("-1: %v\n", A.GetIndex(-1))
    // this should fail: index out of bounds
    //t.Logf(" 6: %v\n", A.GetIndex(6))
    t.Logf(" every 2nd: %v\n", A.GetIndexes(Indexes(0, A.NumElements(), 2)...))
}

func testScalars(t *testing.T) {
    f := FScalar(2.0)
    t.Logf(" f = %v\n", f)
    t.Logf("-f = %v\n", -f)
    z := FScalar(f * f)
    t.Logf(" z = %v\n", z)

}

func TestArrayCreate(t *testing.T) {
    m := FloatVector([]float64{0, 1, 2, 3, 4, 5, 6, 7, 8, 9})
    b := FloatVector(m.FloatArray()[2:5])
    t.Logf("len(m) = %d, len(b) = %d, b=\n%v\n", m.NumElements(), b.NumElements(), b)
}

func TestParseSpe(t *testing.T) {
    s := "{2 3 [3.666666666666667, 3.142857142857143, 4.857142857142857, 4.000000000000000, 5.000000000000000, 6.000000000000000]}"
    m, err := FloatParse(s)
    if err != nil {
        t.Logf("parse error: %v\n", err)
    } else {
        t.Logf("rows: %d, cols: %d, data:\n%v\n", m.Rows(), m.Cols(), m)
    }

    s2 := "{0 1 []}"
    m, err = FloatParse(s2)
    if err != nil {
        t.Logf("parse error: %v\n", err)
    } else {
        t.Logf("rows: %d, cols: %d, data:\n%v\n", m.Rows(), m.Cols(), m)
    }
}

func TestStacked(t *testing.T) {
    a := FloatZeros(3, 3)
    b := FloatOnes(3, 3)
    m0, _ := FloatMatrixStacked(StackDown, a, b)
    m1, _ := FloatMatrixStacked(StackRight, a, b)
    t.Logf("stack down=\n%v\n", m0.ToString("%.2f"))
    t.Logf("stack right=\n%v\n", m1.ToString("%.2f"))
}

func TestFromTable(t *testing.T) {
    data := [][]float64{
        []float64{1, 2, 3},
        []float64{4, 5, 6},
        []float64{7, 8, 9}}

    a := FloatMatrixFromTable(data, RowOrder)
    b := FloatMatrixFromTable(data, ColumnOrder)
    t.Logf("a=\n%v\n", a.ToString("%.2f"))
    t.Logf("b=\n%v\n", b.ToString("%.2f"))
    t.Logf("b == a:   %v\n", b.Equal(a))
    t.Logf("b == a.T: %v\n", b.Equal(a.Transpose()))
}

func TestFromTable2(t *testing.T) {
    data := []float64{1, 2, 3, 4, 5, 6}
    a := FloatNew(2, 3, data)
    b := FloatNew(2, 3, data, RowOrder)
    t.Logf("a=\n%v\n", a)
    t.Logf("b=\n%v\n", b)
}

func TestSubMatrix(t *testing.T) {
    data := [][]float64{
        []float64{1, 2, 3},
        []float64{4, 5, 6},
        []float64{7, 8, 9}}

    A := FloatMatrixFromTable(data)
    b := A.SubMatrix(1, 1)
    t.Logf("A=\n%v\n", A)
    t.Logf("b.NumElements: %d, b=\n%v\n", b.NumElements(), b)
    r0 := A.SubMatrix(0, 0, 1, A.Cols())
    c0 := A.SubMatrix(0, 0, A.Rows(), 1)
    t.Logf("r0=\n%v\n", r0)
    t.Logf("c0=\n%v\n", c0)
    b.Scale(2.0)
    t.Logf("2*b=\n%v\n", b)
    r0.Add(5.0)
    c0.Add(3.0, 1, 2)
    t.Logf("r0+5=\n%v\n", r0)
    t.Logf("c0[1,2]+3.0=\n%v\n", c0)

    t.Logf("A=\n%v\n", A)
    b.SubMatrixOf(A, 0, 0, 2, 2)
    t.Logf("b.SubmatrixOf(A, 0, 0, 2, 2)\n%v\n", b)
    c := b.Copy()
    t.Logf("c = b.Copy()\n%v\n", c)
    c.Exp()
    t.Logf("c.Exp()\n%v\n", c)
    b.Set(c)
    t.Logf("b.Set(c)\n%v\n", b)
    t.Logf("A=\n%v\n", A)
    diag := A.SubMatrix(0, 0, 1, A.Cols(), A.LeadingIndex()+1)
    t.Logf("diag= %d elements\n%v\n", diag.NumElements(), diag)
    diag.SetIndexesFromArray([]float64{1.0, 2.0, 3.0}, 0, 1, 2)
    t.Logf("diag=\n%v\n", diag)
    t.Logf("A=\n%v\n", A)
    foo := A.SubMatrix(0, 0, A.Rows()-1, 2, A.LeadingIndex()+1)
    t.Logf("foo=\n%v\n", foo)
    faa := A.SubMatrix(0, 1, A.Rows()-1, 2, A.LeadingIndex()+1)
    t.Logf("faa=\n%v\n", faa)
    d2 := A.Diag()
    t.Logf("A.Diag()=\n%v\n", d2)
}

// Local Variables:
// tab-width: 4
// End:
