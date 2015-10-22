package symexpr

import (
  "testing"
  
  "fmt"
)


func Test_Print(t *testing.T) {
  fmt.Printf( "Testing: Print                " )
  
  c := NewConstantF(1.57)
  s := NewSin(c)
  t.Logf( "Hello Me!   %v\n", s )

  c.F = 90.0
  t.Logf( "Hello Me!   %v\n", c )
  
  fmt.Printf( "Passed\n" )
}



// func TEST() {
//   A1 := SE.NewAdd()
//   A2 := SE.NewAdd()
//
//   M1 := SE.NewMul()
//   M2 := SE.NewMul()
//   M3 := SE.NewMul()
//
//   c1 := new(SE.Constant)
//   c1.P = 0
//   c2 := new(SE.Constant)
//   c2.P = 1
//   c3 := new(SE.Constant)
//   c3.P = 1
//
//   x1 := new(SE.Var)
//   x1.P = 0
//   x2 := new(SE.Var)
//   x2.P = 0
//
//   M1.Insert(c1)
//   M1.Insert(x1)
//   M2.Insert(c2)
//   M2.Insert(x2)
//   M3.Insert(c3)
//   M3.Insert(x1)
//
//   A1.Insert(M1)
//   A1.Insert(M2)
//
//   A2.Insert(M2)
//   A2.Insert(M3)
//
//
//   m12 := M1.AmIAlmostSame(M2)
//   m13 := M1.AmIAlmostSame(M3)
//   m23 := M2.AmIAlmostSame(M3)
//
//   fmt.Printf( "m12: %v\n", m12 )
//   fmt.Printf( "m12: %v\n", m13 )
//   fmt.Printf( "m12: %v\n", m23 )
//   fmt.Printf( "\n" )
//
//   fmt.Printf( "A1:  %v\n", A1 )
//   fmt.Printf( "A2:  %v\n", A2 )
//   A1.Sort()
//   A2.Sort()
//   fmt.Printf( "A1:  %v\n", A1 )
//   fmt.Printf( "A2:  %v\n", A2 )
// }
/*

func TestAmIAlmostSame() {
	var add1,add2,add3,add4 Add
	add1.CS[0] = &ConstantF{F: 4.20}
	add1.CS[1] = &Var{P:1}
	add2.CS[0] = &ConstantF{F: 4.20}
	add2.CS[1] = &Var{P:1}
	add3.CS[0] = &ConstantF{F: 4.20}
	add3.CS[1] = &Var{P:2}
	add4.CS[0] = &ConstantF{F: -4.20}
	add4.CS[1] = &Var{P:1}

	var a1,a2,a3,a4 Expr
	a1,a2,a3,a4 = &add1,&add2,&add3,&add4

	fmt.Printf( "Testing Almost\n\n" )

	fmt.Printf( "a1: %v\na2: %v\na3: %v\na4: %v\n", a1,a2,a3,a4 )

	fmt.Printf( "a1->a2 = %v (true)\n", a1.AmIAlmostSame(a2) )
	fmt.Printf( "a1->a3 = %v (false)\n", a1.AmIAlmostSame(a3) )
	fmt.Printf( "a1->a4 = %v (true)\n", a1.AmIAlmostSame(a4) )

}

func testAddSimp(srules SimpRules) {
  var add,bdd,cdd Add
  var mul Mul
  mul.CS[0] = &ConstantF{F: 4.20}
  mul.CS[1] = &Var{P:1}
  add.CS[0] = &ConstantF{F: 4.20}
  add.CS[1] = &cdd
  add.CS[2] = &bdd
  add.CS[3] = &mul
  bdd.CS[0] = &Var{P:1}
  bdd.CS[1] = &Var{P:0}
  cdd.CS[0] = &Var{P:1}
  cdd.CS[1] = &Var{P:1}

  var e Expr
  e = &add
//   f = &cdd

  fmt.Printf( "Orig:  %v\n", e )
  es := e.Simplify( srules )
  fmt.Printf( "Simp:  %v\n\n\n", es )
  es.CalcExprStats(0)

/*  for i:= 0; i < es.Size(); i++ {
    p,q := i,i
    fmt.Printf( "%v: %v\n", i, es.GetExpr(&p) )
    et := es.Clone()
    fmt.Printf( "Orig:  %v\n", et )
    f := cdd.Clone()
    eu := SwapExpr( et, f, q )
    if eu != nil {
      et = eu
    } else {
      fmt.Printf( "nil\n" )
    }
    fmt.Printf( "Swap:  %v\n", et )
    ev := et.Simplify(srules)
    fmt.Printf( "Simp:  %v\n\n", ev )
  }*/



//}*/


