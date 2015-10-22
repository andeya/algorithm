package symexpr

import (
	"testing"

	"fmt"
	"math"
)

func Test_SimpLeafs(TEST *testing.T) {
	fmt.Printf("Testing: Simplifying Leafs    ")
	var rules = DefaultRules()

	t := NewTime()
	t_simp := t.Simplify(rules)
	if t != t_simp {
		TEST.Fatalf("FAIL Time: t != t_simp  ~  %v -> %v", t, t_simp)
	} else {
		TEST.Logf("Time Simp:  %v -> %v", t, t_simp)
	}

	v := NewVar(0)
	v_simp := v.Simplify(rules)
	if v != v_simp {
		TEST.Fatalf("FAIL Var: v != v_simp ~  %v ->%v", v, v_simp)
	} else {
		TEST.Logf("Var Simp:  %v -> %v", v, v_simp)
	}

	c := NewConstant(0)
	c_simp := c.Simplify(rules)
	if c != c_simp {
		TEST.Fatalf("FAIL Constant: c != c_simp  ~  %v -> %v", c, c_simp)
	} else {
		TEST.Logf("Constant Simp  %v -> %v", c, c_simp)
	}

	s := NewSystem(0)
	s_simp := s.Simplify(rules)
	if s != s_simp {
		TEST.Fatalf("FAIL System: s != s_simp  ~  %v -> %v", s, s_simp)
	} else {
		TEST.Logf("System Simp:  %v -> %v", s, s_simp)
	}

	f := NewConstantF(3.141529)
	f_simp := f.Simplify(rules)
	if f != f_simp {
		TEST.Fatalf("FAIL ConstantF: f != f_simp  ~  %v -> %v", f, f_simp)
	} else {
		TEST.Logf("ConstantF Simp:  %v -> %v", f, f_simp)
	}

	n := NewConstantF(math.NaN())
	n_simp := n.Simplify(rules)
	if n_simp != nil {
		TEST.Fatalf("FAIL NaN: n_simp != nil  ~  %v -> %v", n, n_simp)
	} else {
		TEST.Logf("NaN Simp:  %v -> %v", n, n_simp)
	}

	i := NewConstantF(math.Inf(1))
	i_simp := i.Simplify(rules)
	if i_simp != nil {
		TEST.Fatalf("FAIL Inf: i == i_simp  ~  %v -> %v", i, i_simp)
	} else {
		TEST.Logf("Inf Simp:  %v -> %v", i, i_simp)
	}

	fmt.Printf("Passed\n")
}

func Test_SimpNeg(TEST *testing.T) {
	fmt.Printf("Testing: Simplifying Neg      ")
	var rules = DefaultRules()

	nn := NewNeg(nil)
	nn_simp := nn.Simplify(rules)
	if nn_simp != nil {
		TEST.Fatalf("FAIL NegNil: nn_simp != nil  ~  %v -> %v", nn, nn_simp)
	} else {
		TEST.Logf("NegNil:  %v -> %v", nn, nn_simp)
	}

	nl := NewNeg(NewNull())
	nl_simp := nl.Clone().Simplify(rules)
	if nl_simp != nil {
		TEST.Fatalf("FAIL NegNull: nn_simp != Null  ~  %v -> %v", nl, nl_simp)
	} else {
		TEST.Logf("NegNull:  %v -> %v", nl, nl_simp)
	}

	nv := NewNeg(NewVar(0))
	nv_simp := nv.Simplify(rules)
	if nv != nv_simp {
		TEST.Fatalf("FAIL NegVar: nv != nv_simp  ~  %v -> %v", nv, nv_simp)
	} else {
		TEST.Logf("NegVar:  %v -> %v", nv, nv_simp)
	}

	c := NewConstant(0)
	nc := NewNeg(c)
	nc_simp := nc.Clone().Simplify(rules)
	if !c.AmISame(nc_simp) {
		TEST.Fatalf("FAIL NegConst: c != nc_simp  ~  %v -> %v", nc, nc_simp)
	} else {
		TEST.Logf("NegConst:  %v -> %v", nc, nc_simp)
	}

	m := NewMul()
	m.Insert(NewVar(0))
	m.Insert(NewConstant(0))
	nm := NewNeg(m)
	nm_simp := nm.Clone().Simplify(rules)
	if !nm_simp.AmISame(m) {
		TEST.Fatalf("FAIL NegMul: m != nm_simp  ~  %v -> %v", nm, nm_simp)
	} else {
		TEST.Logf("NegMul:  %v -> %v", nm, nm_simp)
	}

	fmt.Printf("Passed\n")
}

func Test_SimpAbs(TEST *testing.T) {
	fmt.Printf("Testing: Simplifying Abs      ")
	var rules = DefaultRules()

	nn := NewAbs(nil)
	nn_simp := nn.Clone().Simplify(rules)
	if nn_simp != nil {
		TEST.Fatalf("FAIL AbsNil: nn_simp != nil  ~  %v -> %v", nn, nn_simp)
	} else {
		TEST.Logf("AbsNil:  %v -> %v", nn, nn_simp)
	}

	nl := NewAbs(NewNull())
	nl_simp := nl.Clone().Simplify(rules)
	if nl_simp != nil {
		TEST.Fatalf("FAIL AbsNull: nn_simp != Null  ~  %v -> %v", nl, nl_simp)
	} else {
		TEST.Logf("AbsNull:  %v -> %v", nl, nl_simp)
	}

	nv := NewAbs(NewVar(0))
	nv_simp := nv.Simplify(rules)
	if !nv.AmISame(nv_simp) {
		TEST.Fatalf("FAIL AbsVar: nv != nv_simp  ~  %v -> %v", nv, nv_simp)
	} else {
		TEST.Logf("AbsVar:  %v -> %v", nv, nv_simp)
	}

	c := NewConstant(0)
	nc := NewAbs(c)
	nc_simp := nc.Simplify(rules)
	if !nc.AmISame(nc_simp) {
		TEST.Fatalf("FAIL AbsConst: c != nc_simp  ~  %v -> %v", nc, nc_simp)
	} else {
		TEST.Logf("AbsConst:  %v -> %v", nc, nc_simp)
	}

	m := NewMul()
	m.Insert(NewVar(0))
	m.Insert(NewConstant(0))
	nm := NewAbs(m)
	nm_simp := nm.Simplify(rules)
	if !nm_simp.AmISame(nm) {
		TEST.Fatalf("FAIL AbsMul: m != nm_simp  ~  %v -> %v", nm, nm_simp)
	} else {
		TEST.Logf("AbsMul:  %v -> %v", nm, nm_simp)
	}

	fmt.Printf("Passed\n")
}

func Test_SimpSqrt(TEST *testing.T) {
	fmt.Printf("Testing: Simplifying Sqrt     ")
	var rules = DefaultRules()

	nn := NewSqrt(nil)
	nn_simp := nn.Clone().Simplify(rules)
	if nn_simp != nil {
		TEST.Fatalf("FAIL SqrtNil: nn_simp != nil  ~  %v -> %v", nn, nn_simp)
	} else {
		TEST.Logf("SqrtNil:  %v -> %v", nn, nn_simp)
	}

	nl := NewSqrt(NewNull())
	nl_simp := nl.Clone().Simplify(rules)
	if nl_simp != nil {
		TEST.Fatalf("FAIL SqrtNull: nn_simp != Null  ~  %v -> %v", nl, nl_simp)
	} else {
		TEST.Logf("SqrtNull:  %v -> %v", nl, nl_simp)
	}

	ns := NewSqrt(NewSqrt(NewVar(0)))
	ns_simp := ns.Clone().Simplify(rules)
	ns_corr := NewPowF(NewVar(0), 0.25)
	if !ns_simp.AmISame(ns_corr) {
		TEST.Fatalf("FAIL SqrtSqrt: ns_simp != ns_corr  ~  %v -> %v  ==  %v", ns, ns_simp, ns_corr)
	} else {
		TEST.Logf("SqrtSqrt:  %v -> %v", ns, ns_simp)
	}

	c := NewConstant(0)
	nc := NewSqrt(c)
	nc_simp := nc.Clone().Simplify(rules)
	if !c.AmISame(nc_simp) {
		TEST.Fatalf("FAIL SqrtConst: c != nc_simp  ~  %v -> %v", nc, nc_simp)
	} else {
		TEST.Logf("SqrtConst:  %v -> %v", nc, nc_simp)
	}

	nf := NewSqrt(NewPowF(NewVar(0), 3.14))
	nf_simp := nf.Clone().Simplify(rules)
	nf_corr := NewPowF(NewVar(0), 3.14/2.0)
	if !nf_simp.AmISame(nf_corr) {
		TEST.Fatalf("FAIL SqrtPowF: nf_simp != nf_corr  ~  %v -> %v  ==  %v", nf, nf_simp, nf_corr)
	} else {
		TEST.Logf("SqrtPowF:  %v -> %v", nf, nf_simp)
	}

	ni1 := NewSqrt(NewPowI(NewVar(0), 1))
	ni1_simp := ni1.Clone().Simplify(rules)
	ni1_corr := NewSqrt(NewVar(0))
	if !ni1_simp.AmISame(ni1_corr) {
		TEST.Fatalf("FAIL SqrtPowI(ni1): ni1_simp != ni1_corr  ~  %v -> %v  ==  %v", ni1, ni1_simp, ni1_corr)
	} else {
		TEST.Logf("SqrtPowI(ni1):  %v -> %v", ni1, ni1_simp)
	}

	//   ni1n := NewSqrt(NewPowI(NewVar(0),-1))
	//   ni1n_simp := ni1n.Clone().Simplify(rules)
	//   ni1n_corr := NewSqrt(NewVar(0))
	//   if !ni1n_simp.AmISame(ni1n_corr) {
	// 	TEST.Fatalf( "FAIL SqrtPowI(ni1n): ni1n_simp != ni1n_corr  ~  %v -> %v  ==  %v", ni1n,ni1n_simp, ni1n_corr )
	//   } else {
	// 	TEST.Logf( "SqrtPowI(ni1n):  %v -> %v", ni1n,ni1n_simp )
	//   }
	//
	//   ni2 := NewSqrt(NewPowI(NewVar(0),2))
	//   ni2_simp := ni2.Clone().Simplify(rules)
	//   ni2_corr := NewSqrt(NewVar(0))
	//   if !ni2_simp.AmISame(ni2_corr) {
	// 	TEST.Fatalf( "FAIL SqrtPowI(ni2): ni2_simp != ni2_corr  ~  %v -> %v  ==  %v", ni2,ni2_simp, ni2_corr )
	//   } else {
	// 	TEST.Logf( "SqrtPowI(ni2):  %v -> %v", ni2,ni2_simp )
	//   }
	//
	//   ni2n := NewSqrt(NewPowI(NewVar(0),-2))
	//   ni2n_simp := ni2n.Clone().Simplify(rules)
	//   ni2n_corr := NewSqrt(NewVar(0))
	//   if !ni2n_simp.AmISame(ni2n_corr) {
	// 	TEST.Fatalf( "FAIL SqrtPowI(ni2n): ni2n_simp != ni2n_corr  ~  %v -> %v  ==  %v", ni2n,ni2n_simp, ni2n_corr )
	//   } else {
	// 	TEST.Logf( "SqrtPowI(ni2n):  %v -> %v", ni2n,ni2n_simp )
	//   }
	//
	//   niOdd := NewSqrt(NewPowI(NewVar(0),3))
	//   niOdd_simp := niOdd.Clone().Simplify(rules)
	//   niOdd_corr := NewSqrt(NewVar(0))
	//   if !niOdd_simp.AmISame(niOdd_corr) {
	// 	TEST.Fatalf( "FAIL SqrtPowI(niOdd): niOdd_simp != niOdd_corr  ~  %v -> %v  ==  %v", niOdd,niOdd_simp, niOdd_corr )
	//   } else {
	// 	TEST.Logf( "SqrtPowI(niOdd):  %v -> %v", niOdd,niOdd_simp )
	//   }
	//
	//   niOddn := NewSqrt(NewPowI(NewVar(0),-3))
	//   niOddn_simp := niOddn.Clone().Simplify(rules)
	//   niOddn_corr := NewSqrt(NewVar(0))
	//   if !niOddn_simp.AmISame(niOddn_corr) {
	// 	TEST.Fatalf( "FAIL SqrtPowI(niOddn): niOddn_simp != niOddn_corr  ~  %v -> %v  ==  %v", niOddn,niOddn_simp, niOddn_corr )
	//   } else {
	// 	TEST.Logf( "SqrtPowI(niOddn):  %v -> %v", niOddn,niOddn_simp )
	//   }
	//
	//   ni2N := NewSqrt(NewPowI(NewVar(0),4))
	//   ni2N_simp := ni2N.Clone().Simplify(rules)
	//   ni2N_corr := NewSqrt(NewVar(0))
	//   if !ni2N_simp.AmISame(ni2N_corr) {
	// 	TEST.Fatalf( "FAIL SqrtPowI(ni2N): ni2N_simp != ni2N_corr  ~  %v -> %v  ==  %v", ni2N,ni2N_simp, ni2N_corr )
	//   } else {
	// 	TEST.Logf( "SqrtPowI(ni2N):  %v -> %v", ni2N,ni2N_simp )
	//   }
	//
	//   ni2Nn := NewSqrt(NewPowI(NewVar(0),-4))
	//   ni2Nn_simp := ni2Nn.Clone().Simplify(rules)
	//   ni2Nn_corr := NewSqrt(NewVar(0))
	//   if !ni2Nn_simp.AmISame(ni2Nn_corr) {
	// 	TEST.Fatalf( "FAIL SqrtPowI(ni2Nn): ni2Nn_simp != ni2Nn_corr  ~  %v -> %v  ==  %v", ni2Nn,ni2Nn_simp, ni2Nn_corr )
	//   } else {
	// 	TEST.Logf( "SqrtPowI(ni2Nn):  %v -> %v", ni2Nn,ni2Nn_simp )
	//   }
	//
	//

	// the following should not change

	nv := NewSqrt(NewVar(0))
	nv_simp := nv.Simplify(rules)
	if nv != nv_simp {
		TEST.Fatalf("FAIL SqrtVar: nv != nv_simp  ~  %v -> %v", nv, nv_simp)
	} else {
		TEST.Logf("SqrtVar:  %v -> %v", nv, nv_simp)
	}

	m := NewMul()
	m.Insert(NewVar(0))
	m.Insert(NewConstant(0))
	nm := NewSqrt(m)
	nm_simp := nm.Simplify(rules)
	if nm_simp != nm {
		TEST.Fatalf("FAIL SqrtMul: nm != nm_simp  ~  %v -> %v", nm, nm_simp)
	} else {
		TEST.Logf("SqrtMul:  %v -> %v", nm, nm_simp)
	}

	fmt.Printf("Passed\n")
}

func Test_SimpAdd(TEST *testing.T) {
	fmt.Printf("Testing: Simplifying Add     ")
	var rules = DefaultRules()

	an := NewAdd()
	an.Insert(nil)
	an_simp := an.Clone().Simplify(rules)
	if an_simp != nil {
		TEST.Fatalf("FAIL AddNil: an_simp != nil  ~  %v -> %v", an, an_simp)
	} else {
		TEST.Logf("AddNil:  %v -> %v", an, an_simp)
	}

	al := NewAdd()
	al.Insert(NewNull())
	al_simp := al.Clone().Simplify(rules)
	if al_simp != nil {
		TEST.Fatalf("FAIL AddNull: nn_simp != Null  ~  %v -> %v", al, al_simp)
	} else {
		TEST.Logf("AddNull:  %v -> %v", al, al_simp)
	}

	as := NewAdd()
	as.Insert(NewVar(0))
	as_simp := as.Clone().Simplify(rules)
	as_corr := NewVar(0)
	if !as_simp.AmISame(as_corr) {
		TEST.Fatalf("FAIL AddVar: as_simp != as_corr  ~  %v -> %v  ==  %v", as, as_simp, as_corr)
	} else {
		TEST.Logf("AddVar:  %v -> %v", as, as_simp)
	}

	c := NewConstant(0)
	ac := NewAdd()
	ac.Insert(c)
	ac_simp := ac.Clone().Simplify(rules)
	if !c.AmISame(ac_simp) {
		TEST.Fatalf("FAIL AddConst: c != ac_simp  ~  %v -> %v", ac, ac_simp)
	} else {
		TEST.Logf("AddConst:  %v -> %v", ac, ac_simp)
	}

	// x + y
	A1 := NewAdd()
	A1.Insert(NewVar(0))
	A1.Insert(NewVar(1))
	// A1.CS[0], A1.CS[1] = A1.CS[1], A1.CS[0]
	A1_simp := A1.Clone().Simplify(rules)
	A1_corr := A1.Clone()
	if !A1_simp.AmISame(A1_corr) {
		TEST.Fatalf("FAIL Add_A1: A1_simp != A1_corr  ~  %v -> %v  ==  %v", A1, A1_simp, A1_corr)
	} else {
		TEST.Logf("AddVar:  %v -> %v", A1, A1_simp)
	}

	// x + x
	A2 := NewAdd()
	A2.Insert(NewVar(0))
	A2.Insert(NewVar(0))
	A2_simp := A2.Clone().Simplify(rules)
	A2_corr := NewMul()
	A2_corr.Insert(NewConstant(0))
	A2_corr.Insert(NewVar(0))

	if !A2_simp.AmISame(A2_corr) {
		TEST.Fatalf("FAIL Add_A2: A2_simp != A2_corr  ~  %v -> %v  ==  %v", A2, A2_simp, A2_corr)
	} else {
		TEST.Logf("AddVar:  %v -> %v", A2, A2_simp)
	}

	// x + x + y
	A3 := NewAdd()
	A3.Insert(NewVar(0))
	A3.Insert(NewVar(0))
	A3.Insert(NewVar(1))
	A3_simp := A3.Clone().Simplify(rules)
	A3_corr := NewAdd()
	A3m := NewMul()
	A3m.Insert(NewConstant(0))
	A3m.Insert(NewVar(0))
	A3_corr.Insert(A3m)
	A3_corr.Insert(NewVar(1))
	if !A3_simp.AmISame(A3_corr) {
		TEST.Fatalf("FAIL Add_A3: A3_simp != A3_corr  ~  %v -> %v  ==  %v", A3, A3_simp, A3_corr)
	} else {
		TEST.Logf("AddVar:  %v -> %v", A3, A3_simp)
	}

	// x + ( x + y )
	A4 := NewAdd()
	A4.Insert(NewVar(0))
	A4b := NewAdd()
	A4b.Insert(NewVar(0))
	A4b.Insert(NewVar(1))
	A4.Insert(A4b)
	A4_simp := A4.Clone().Simplify(rules)
	A4_corr := NewAdd()
	A4m := NewMul()
	A4m.Insert(NewConstant(0))
	A4m.Insert(NewVar(0))
	A4_corr.Insert(A4m)
	A4_corr.Insert(NewVar(1))
	if !A4_simp.AmISame(A4_corr) {
		TEST.Fatalf("FAIL Add_A4: A4_simp != A4_corr  ~  %v -> %v  ==  %v", A4, A4_simp, A4_corr)
	} else {
		TEST.Logf("AddVar:  %v -> %v", A4, A4_simp)
	}

	//  F(x) +  F(x)
	// x^2 + x^2
	A5 := NewAdd()
	A5.Insert(NewPowI(NewVar(0), 2))
	A5.Insert(NewPowI(NewVar(0), 2))
	A5_simp := A5.Clone().Simplify(rules)
	A5_corr := NewMul()
	A5_corr.Insert(NewConstant(0))
	A5_corr.Insert(NewPowI(NewVar(0), 2))
	if !A5_simp.AmISame(A5_corr) {
		TEST.Fatalf("FAIL Add_A5: A5_simp != A5_corr  ~  %v -> %v  ==  %v", A5, A5_simp, A5_corr)
	} else {
		TEST.Logf("AddVar:  %v -> %v", A5, A5_simp)
	}

	// cF(x) + dF(x)
	// cx^3.14 + dx^3.15
	A6 := NewAdd()
	A6m1 := NewMul()
	A6m2 := NewMul()
	A6m1.Insert(NewConstant(0))
	A6m1.Insert(NewPowF(NewVar(0), 3.14))
	A6m2.Insert(NewConstant(0))
	A6m2.Insert(NewPowF(NewVar(0), 3.14))
	A6.Insert(A6m1)
	A6.Insert(A6m2)
	A6_simp := A6.Clone().Simplify(rules)
	A6_corr := NewMul()
	A6_corr.Insert(NewConstant(0))
	A6_corr.Insert(NewPowF(NewVar(0), 3.14))
	if !A6_simp.AmISame(A6_corr) {
		TEST.Fatalf("FAIL Add_A6: A6_simp != A6_corr  ~  %v -> %v  ==  %v", A6, A6_simp, A6_corr)
	} else {
		TEST.Logf("AddVar:  %v -> %v", A6, A6_simp)
	}

	// cF(x) +  F(x)
	// cExp(x) + Exp(x)
	A7 := NewAdd()
	A7m1 := NewMul()
	A7m1.Insert(NewConstant(0))
	A7m1.Insert(NewExp(NewVar(0)))
	A7.Insert(A7m1)
	A7.Insert(NewExp(NewVar(0)))
	A7_simp := A7.Clone().Simplify(rules)
	A7_corr := NewMul()
	A7_corr.Insert(NewConstant(0))
	A7_corr.Insert(NewExp(NewVar(0)))
	if !A7_simp.AmISame(A7_corr) {
		TEST.Fatalf("FAIL Add_A7: A7_simp != A7_corr  ~  %v -> %v  ==  %v", A7, A7_simp, A7_corr)
	} else {
		TEST.Logf("AddVar:  %v -> %v", A7, A7_simp)
	}

	//  F(x) + cF(x)
	//  sin(x) + csin(x)
	A8 := NewAdd()
	A8m2 := NewMul()
	A8m2.Insert(NewConstant(0))
	A8m2.Insert(NewSin(NewVar(0)))
	A8.Insert(NewSin(NewVar(0)))
	A8.Insert(A8m2)
	A8_simp := A8.Clone().Simplify(rules)
	A8_corr := NewMul()
	A8_corr.Insert(NewConstant(0))
	A8_corr.Insert(NewSin(NewVar(0)))
	if !A8_simp.AmISame(A8_corr) {
		TEST.Fatalf("FAIL Add_A8: A8_simp != A8_corr  ~  %v -> %v  ==  %v", A8, A8_simp, A8_corr)
	} else {
		TEST.Logf("AddVar:  %v -> %v", A8, A8_simp)
	}

	//  F(x) -  F(x)
	// x - x
	A9 := NewAdd()
	A9.Insert(NewVar(0))
	A9.Insert(NewNeg(NewVar(0)))
	A9_simp := A9.Clone().Simplify(rules)
	A9_corr := NewMul()
	A9_corr.Insert(NewConstant(0))
	A9_corr.Insert(NewVar(0))

	if !A9_simp.AmISame(A9_corr) {
		TEST.Fatalf("FAIL Add_A9: A9_simp != A9_corr  ~  %v -> %v  ==  %v", A9, A9_simp, A9_corr)
	} else {
		TEST.Logf("AddVar:  %v -> %v", A9, A9_simp)
	}
	// x + (-x+y)
	A10 := NewAdd()
	A10.Insert(NewVar(0))
	A10b := NewAdd()
	A10b.Insert(NewNeg(NewVar(0)))
	A10b.Insert(NewVar(1))
	A10.Insert(A10b)
	A10_simp := A10.Clone().Simplify(rules)
	A10_corr := NewAdd()
	A10m := NewMul()
	A10m.Insert(NewConstant(0))
	A10m.Insert(NewVar(0))
	A10_corr.Insert(A10m)
	A10_corr.Insert(NewVar(1))
	if !A10_simp.AmISame(A10_corr) {
		TEST.Fatalf("FAIL Add_A10: A10_simp != A10_corr  ~  %v -> %v  ==  %v", A10, A10_simp, A10_corr)
	} else {
		TEST.Logf("AddVar:  %v -> %v", A10, A10_simp)
	}
	// x^2 - x^2
	A11 := NewAdd()
	A11.Insert(NewPowI(NewVar(0), 2))
	A11.Insert(NewNeg(NewPowI(NewVar(0), 2)))
	A11_simp := A11.Clone().Simplify(rules)
	A11_corr := NewMul()
	A11_corr.Insert(NewConstant(0))
	A11_corr.Insert(NewPowI(NewVar(0), 2))
	if !A11_simp.AmISame(A11_corr) {
		TEST.Fatalf("FAIL Add_A11: A11_simp != A11_corr  ~  %v -> %v  ==  %v", A11, A11_simp, A11_corr)
	} else {
		TEST.Logf("AddVar:  %v -> %v", A11, A11_simp)
	}

	// cF(x) - dF(x)
	// cx^3.14 - dx^3.15
	A12 := NewAdd()
	A12m1 := NewMul()
	A12m2 := NewMul()
	A12m1.Insert(NewConstant(0))
	A12m1.Insert(NewPowF(NewVar(0), 3.14))
	A12m2.Insert(NewConstant(0))
	A12m2.Insert(NewPowF(NewVar(0), 3.14))
	A12.Insert(A12m1)
	A12.Insert(NewNeg(A12m2))
	A12_simp := A12.Clone().Simplify(rules)
	A12_corr := NewMul()
	A12_corr.Insert(NewConstant(0))
	A12_corr.Insert(NewPowF(NewVar(0), 3.14))
	if !A12_simp.AmISame(A12_corr) {
		TEST.Fatalf("FAIL Add_A12: A12_simp != A12_corr  ~  %v -> %v  ==  %v", A12, A12_simp, A12_corr)
	} else {
		TEST.Logf("AddVar:  %v -> %v", A12, A12_simp)
	}

	// cF(x) -  F(x)
	// cExp(x) - Exp(x)
	A13 := NewAdd()
	A13m1 := NewMul()
	A13m1.Insert(NewConstant(0))
	A13m1.Insert(NewExp(NewVar(0)))
	A13.Insert(A13m1)
	A13.Insert(NewNeg(NewExp(NewVar(0))))
	A13_simp := A13.Clone().Simplify(rules)
	A13_corr := NewMul()
	A13_corr.Insert(NewConstant(0))
	A13_corr.Insert(NewExp(NewVar(0)))
	if !A13_simp.AmISame(A13_corr) {
		TEST.Fatalf("FAIL Add_A13: A13_simp != A13_corr  ~  %v -> %v  ==  %v", A13, A13_simp, A13_corr)
	} else {
		TEST.Logf("AddVar:  %v -> %v", A13, A13_simp)
	}

	//  F(x) - cF(x)
	//  sin(x) - csin(x)
	A14 := NewAdd()
	A14m2 := NewMul()
	A14m2.Insert(NewConstant(0))
	A14m2.Insert(NewSin(NewVar(0)))
	A14.Insert(NewSin(NewVar(0)))
	A14.Insert(NewNeg(A14m2))
	A14_simp := A14.Clone().Simplify(rules)
	A14_corr := NewMul()
	A14_corr.Insert(NewConstant(0))
	A14_corr.Insert(NewSin(NewVar(0)))
	if !A14_simp.AmISame(A14_corr) {
		TEST.Fatalf("FAIL Add_A14: A14_simp != A14_corr  ~  %v -> %v  ==  %v", A14, A14_simp, A14_corr)
	} else {
		TEST.Logf("AddVar:  %v -> %v", A14, A14_simp)
	}

	// c0 + c1*x + c2*x*x*x
	A15 := NewAdd()
	A15.Insert(NewConstant(0))
	A15m1 := NewMul()
	A15m1.Insert(NewConstant(1))
	A15m1.Insert(NewVar(0))
	A15.Insert(A15m1)
	A15m2 := NewMul()
	A15m2.Insert(NewConstant(2))
	A15m2.Insert(NewVar(0))
	A15m2.Insert(NewVar(0))
	A15m2.Insert(NewVar(0))
	A15.Insert(A15m2)
	A15_simp := A15.Clone().Simplify(rules)
	A15_corr := NewAdd()
	A15_corr.Insert(NewConstant(0))
	A15_corrm1 := NewMul()
	A15_corrm1.Insert(NewConstant(1))
	A15_corrm1.Insert(NewVar(0))
	A15_corr.Insert(A15_corrm1)
	A15_corrm2 := NewMul()
	A15_corrm2.Insert(NewConstant(2))
	A15_corrm2.Insert(NewPowI(NewVar(0), 3))
	A15_corr.Insert(A15_corrm2)
	if !A15_simp.AmISame(A15_corr) {
		TEST.Fatalf("FAIL Add_A15: A15_simp != A15_corr  ~  %v -> %v  ==  %v", A15, A15_simp, A15_corr)
	} else {
		TEST.Logf("AddVar:  %v -> %v", A15, A15_simp)
	}

	// X - (X+Y) ???

}

func Test_SimpMul(TEST *testing.T) {
	fmt.Printf("Testing: Simplifying Mul      ")
	var rules = DefaultRules()

	mn := NewMul()
	mn.Insert(nil)
	mn_simp := mn.Clone().Simplify(rules)
	if mn_simp != nil {
		TEST.Fatalf("FAIL MulNil: mn_simp != nil  ~  %v -> %v", mn, mn_simp)
	} else {
		TEST.Logf("MulNil:  %v -> %v", mn, mn_simp)
	}

	ml := NewMul()
	ml.Insert(NewNull())
	ml_simp := ml.Clone().Simplify(rules)
	if ml_simp != nil {
		TEST.Fatalf("FAIL MulNull: nn_simp != Null  ~  %v -> %v", ml, ml_simp)
	} else {
		TEST.Logf("MulNull:  %v -> %v", ml, ml_simp)
	}

	mv := NewMul()
	mv.Insert(NewVar(0))
	mv_simp := mv.Clone().Simplify(rules)
	mv_corr := NewVar(0)
	if !mv_simp.AmISame(mv_corr) {
		TEST.Fatalf("FAIL MulVar: mv_simp != mv_corr  ~  %v -> %v  ==  %v", mv, mv_simp, mv_corr)
	} else {
		TEST.Logf("MulVar:  %v -> %v", mv, mv_simp)
	}

	c := NewConstant(0)
	mc := NewMul()
	mc.Insert(c)
	mc_simp := mc.Clone().Simplify(rules)
	if !c.AmISame(mc_simp) {
		TEST.Fatalf("FAIL MulConst: c != mc_simp  ~  %v -> %v", mc, mc_simp)
	} else {
		TEST.Logf("MulConst:  %v -> %v", mc, mc_simp)
	}

	cN := NewConstant(0)
	mcN := NewMul()
	mcN.Insert(cN)
	mcN.Insert(NewConstant(1))
	mcN.Insert(NewConstant(2))
	mcN_simp := mcN.Clone().Simplify(rules)
	if !cN.AmISame(mcN_simp) {
		TEST.Fatalf("FAIL MulConst: c != mcN_simp  ~  %v -> %v", mcN, mcN_simp)
	} else {
		TEST.Logf("MulConst:  %v -> %v", mcN, mcN_simp)
	}

	mvN := NewMul()
	mvN.Insert(NewVar(0))
	mvN.Insert(NewVar(0))
	mvN.Insert(NewVar(0))
	mvN_simp := mvN.Clone().Simplify(rules)
	mvN_corr := NewPowI(NewVar(0), 3)
	if !mvN_simp.AmISame(mvN_corr) {
		TEST.Fatalf("FAIL MulVar: mvN_simp != mvN_corr  ~  %v -> %v  ==  %v", mvN, mvN_simp, mvN_corr)
	} else {
		TEST.Logf("MulVar:  %v -> %v", mvN, mvN_simp)
	}

	mcvN := NewMul()
	mcvN.Insert(NewConstant(0))
	mcvN.Insert(NewConstant(1))
	mcvN.Insert(NewConstant(2))
	mcvN.Insert(NewVar(0))
	mcvN.Insert(NewVar(0))
	mcvN.Insert(NewVar(0))
	mcvN_simp := mcvN.Clone().Simplify(rules)
	mcvN_corr := NewMul()
	mcvN_corr.Insert(NewConstant(0))
	mcvN_corr.Insert(NewPowI(NewVar(0), 3))
	if !mcvN_simp.AmISame(mcvN_corr) {
		TEST.Fatalf("FAIL MulVar: mvN_simp != mvN_corr  ~  %v -> %v  ==  %v", mcvN, mcvN_simp, mcvN_corr)
	} else {
		TEST.Logf("MulVar:  %v -> %v", mcvN, mcvN_simp)
	}

	mcvNb := NewMul()
	mcvNb.Insert(NewConstant(0))
	mcvNb.Insert(NewConstant(1))
	mcvNb.Insert(NewVar(0))
	mcvNb.Insert(NewVar(0))
	mcvNb.Insert(NewVar(0))
	mcvNba := NewAdd()
	mcvNba.Insert(NewVar(0))
	mcvNbm := NewMul()
	mcvNbm.Insert(NewConstant(-1))
	mcvNbm.Insert(NewVar(0))
	mcvNba.Insert(mcvNbm)
	mcvNb.Insert(mcvNba)
	mcvNb_simp := mcvNb.Clone().Simplify(rules)
	mcvNb_corr := NewMul()
	mcvNb_corr.Insert(NewConstant(0))
	mcvNb_corr.Insert(NewPowI(NewVar(0), 4))
	if !mcvNb_simp.AmISame(mcvNb_corr) {
		TEST.Fatalf("FAIL MulVar: mvN_simp != mvN_corr  ~  %v -> %v  ==  %v", mcvNb, mcvNb_simp, mcvNb_corr)
	} else {
		TEST.Logf("MulVar:  %v -> %v", mcvNb, mcvNb_simp)
	}

	md := NewMul()
	md.Insert(NewConstant(0))
	mdiv := NewDiv(NewConstant(1), NewVar(0))
	md.Insert(mdiv)
	md_simp := md.Clone().Simplify(rules)
	md_corr := NewDiv(NewConstant(0), NewVar(0))
	if !md_simp.AmISame(md_corr) {
		TEST.Fatalf("FAIL MulVar: md_simp != md_corr  ~  %v -> %v  ==  %v", md, md_simp, md_corr)
	} else {
		TEST.Logf("MulVar:  %v -> %v", md, md_simp)
	}

	md2 := NewMul()
	md2.Insert(NewConstant(0))
	md2m := NewMul()
	md2m.Insert(NewConstant(1))
	md2m.Insert(NewVar(0))
	md2div := NewDiv(md2m, NewVar(0))
	md2.Insert(md2div)
	md2_corr := NewConstant(0)
	md2_simp := md2div.Clone().Simplify(rules)
	if !md2_simp.AmISame(md2_corr) {
		TEST.Fatalf("FAIL MulVar: md2_simp != md2_corr  ~  %v -> %v  ==  %v", md2, md2_simp, md2_corr)
	} else {
		TEST.Logf("MulVar:  %v -> %v", md2, md2_simp)
	}

}

func Test_Damd(TEST *testing.T) {
	fmt.Printf("Testing: Damd Simps      \n")
	var rules = DefaultRules()
	vnames := []string{"x", "y", "z"}

	f1 := parse("2*x^2* x^3", vnames)
	fmt.Printf("f1 orig: %v\n", f1)
	f1s := f1.Clone().Simplify(rules)
	fmt.Printf("f1 simp: %v\n\n", f1s)

	f2 := parse("2*x^-3 * x^2", vnames)
	fmt.Printf("f2 orig: %v\n", f2)
	f2s := f2.Clone().Simplify(rules)
	fmt.Printf("f2 simp: %v\n\n", f2s)

}

func Test_ConstantF(TEST *testing.T) {
	fmt.Printf("Testing: ConstantF Simps\n---------------------\n\n")
	var rules = DefaultRules()
	rules.ConvertConsts = false

	vnames := []string{"x", "y", "z"}

	f1 := parse("x + x", vnames)
	fmt.Printf("f1 orig: %v\n", f1)
	f1s := f1.Clone().Simplify(rules)
	fmt.Printf("f1 simp: %v\n\n", f1s)

	f2 := parse("2*x + x", vnames)
	fmt.Printf("f2 orig: %v\n", f2)
	f2s := f2.Clone().Simplify(rules)
	fmt.Printf("f2 simp: %v\n\n", f2s)

	f3 := parse("x + 2*x", vnames)
	fmt.Printf("f3 orig: %v\n", f3)
	f3s := f3.Clone().Simplify(rules)
	fmt.Printf("f3 simp: %v\n\n", f3s)

	f4 := parse("2*x + 3*x", vnames)
	fmt.Printf("f4 orig: %v\n", f4)
	f4s := f4.Clone().Simplify(rules)
	fmt.Printf("f4 simp: %v\n\n", f4s)

	f5 := parse("x - x", vnames)
	fmt.Printf("f5 orig: %v\n", f5)
	f5s := f5.Clone().Simplify(rules)
	fmt.Printf("f5 simp: %v\n\n", f5s)

	f6 := parse("2*x - x", vnames)
	fmt.Printf("f6 orig: %v\n", f6)
	f6s := f6.Clone().Simplify(rules)
	fmt.Printf("f6 simp: %v\n\n", f6s)

	f7 := parse("x - 2*x", vnames)
	fmt.Printf("f7 orig: %v\n", f7)
	f7s := f7.Clone().Simplify(rules)
	fmt.Printf("f7 simp: %v\n\n", f7s)

	f8 := parse("2*x - 3*x", vnames)
	fmt.Printf("f8 orig: %v\n", f8)
	f8s := f8.Clone().Simplify(rules)
	fmt.Printf("f8 simp: %v\n\n", f8s)

	f9 := parse("5*x - 2*x", vnames)
	fmt.Printf("f9 orig: %v\n", f9)
	f9s := f9.Clone().Simplify(rules)
	fmt.Printf("f9 simp: %v\n\n", f9s)

	/*
		XX := parse("x + x", vnames)
		fmt.Printf("XX orig: %v\n", XX)
		XXs := XX.Clone().Simplify(rules)
		fmt.Printf("XX simp: %v\n\n", XXs)
	*/
}

func Test_Div(TEST *testing.T) {
	fmt.Printf("Testing: Div Simps\n---------------------\n\n")
	var rules = DefaultRules()
	vnames := []string{"x", "y", "z"}
	rules.ConvertConsts = false

	/*****************************************/
	/*****************************************/

	f1 := parse("x / x", vnames)
	f1_corr := NewConstant(-1)
	f1_simp := f1.Clone().Simplify(rules)
	if !f1_simp.AmISame(f1_corr) {
		TEST.Fatalf("FAIL Div: f1_simp != f1_corr  ~  %v -> %v  ==  %v", f1, f1_simp, f1_corr)
	} else {
		fmt.Printf("Div:  %v -> %v", f1, f1_simp)
		TEST.Logf("Div:  %v -> %v", f1, f1_simp)
	}

	f2 := parse("( 2 * x ) / x", vnames)
	f2_corr := NewConstant(-1)
	f2_simp := f2.Clone().Simplify(rules)
	if !f2_simp.AmISame(f2_corr) {
		TEST.Fatalf("FAIL Div: f2_simp != f2_corr  ~  %v -> %v  ==  %v", f2, f2_simp, f2_corr)
	} else {
		fmt.Printf("Div:  %v -> %v", f2, f2_simp)
		TEST.Logf("Div:  %v -> %v", f2, f2_simp)
	}

	f3 := parse("x^2 / x", vnames)
	f3_corr := parse("x", vnames)
	f3_simp := f3.Clone().Simplify(rules)
	if !f3_simp.AmISame(f3_corr) {
		TEST.Fatalf("FAIL Div: f3_simp != f3_corr  ~  %v -> %v  ==  %v", f3, f3_simp, f3_corr)
	} else {
		fmt.Printf("Div:  %v -> %v", f3, f3_simp)
		TEST.Logf("Div:  %v -> %v", f3, f3_simp)
	}

	f4 := parse("x^3 / x^2", vnames)
	f4_corr := parse("x", vnames)
	f4_simp := f4.Clone().Simplify(rules)
	if !f4_simp.AmISame(f4_corr) {
		TEST.Fatalf("FAIL Div: f4_simp != f4_corr  ~  %v -> %v  ==  %v", f4, f4_simp, f4_corr)
	} else {
		fmt.Printf("Div:  %v -> %v", f4, f4_simp)
		TEST.Logf("Div:  %v -> %v", f4, f4_simp)
	}

	f5 := parse("x^4 / x^2", vnames)
	f5_corr := parse("x^2", vnames)
	f5_simp := f5.Clone().Simplify(rules)
	if !f5_simp.AmISame(f5_corr) {
		TEST.Fatalf("FAIL Div: f5_simp != f5_corr  ~  %v -> %v  ==  %v", f5, f5_simp, f5_corr)
	} else {
		fmt.Printf("Div:  %v -> %v", f5, f5_simp)
		TEST.Logf("Div:  %v -> %v", f5, f5_simp)
	}

	f6 := parse("x / x^2", vnames)
	f6_corr := NewDiv(NewConstant(-1), NewVar(0))
	f6_simp := f6.Clone().Simplify(rules)
	if !f6_simp.AmISame(f6_corr) {
		TEST.Fatalf("FAIL Div: f6_simp != f6_corr  ~  %v -> %v  ==  %v", f6, f6_simp, f6_corr)
	} else {
		fmt.Printf("Div:  %v -> %v", f6, f6_simp)
		TEST.Logf("Div:  %v -> %v", f6, f6_simp)
	}

	f7 := parse("x^2 / x^3", vnames)
	f7_corr := NewDiv(NewConstant(-1), NewVar(0))
	f7_simp := f7.Clone().Simplify(rules)
	if !f7_simp.AmISame(f7_corr) {
		TEST.Fatalf("FAIL Div: f7_simp != f7_corr  ~  %v -> %v  ==  %v", f7, f7_simp, f7_corr)
	} else {
		fmt.Printf("Div:  %v -> %v", f7, f7_simp)
		TEST.Logf("Div:  %v -> %v", f7, f7_simp)
	}

	f8 := parse("x^2 / x^4", vnames)
	f8_corr := NewDiv(NewConstant(-1), NewPowI(NewVar(0), 2))
	f8_simp := f8.Clone().Simplify(rules)
	if !f8_simp.AmISame(f8_corr) {
		TEST.Fatalf("FAIL Div: f8_simp != f8_corr  ~  %v -> %v  ==  %v", f8, f8_simp, f8_corr)
	} else {
		fmt.Printf("Div:  %v -> %v", f8, f8_simp)
		TEST.Logf("Div:  %v -> %v", f8, f8_simp)
	}

	/*****************************************/
	/*****************************************/

	f1_n := parse("(y * x) / x", vnames)
	f1_n_corr := NewVar(1)
	f1_n_simp := f1_n.Clone().Simplify(rules)
	if !f1_n_simp.AmISame(f1_n_corr) {
		TEST.Fatalf("FAIL Div: f1_n_simp != f1_n_corr  ~  %v -> %v  ==  %v", f1_n, f1_n_simp, f1_n_corr)
	} else {
		fmt.Printf("Div:  %v -> %v", f1_n, f1_n_simp)
		TEST.Logf("Div:  %v -> %v", f1_n, f1_n_simp)
	}

	f2_n := parse("( 2 * x * y ) / x", vnames)
	f2_n_corr := parse("2 * y", vnames)
	f2_n_simp := f2_n.Clone().Simplify(rules)
	if !f2_n_simp.AmISame(f2_n_corr) {
		TEST.Fatalf("FAIL Div: f2_n_simp != f2_n_corr  ~  %v -> %v  ==  %v", f2_n, f2_n_simp, f2_n_corr)
	} else {
		fmt.Printf("Div:  %v -> %v", f2_n, f2_n_simp)
		TEST.Logf("Div:  %v -> %v", f2_n, f2_n_simp)
	}

	f3_n := parse("( y * x^2 ) / x", vnames)
	f3_n_corr := parse("y * x", vnames)
	f3_n_simp := f3_n.Clone().Simplify(rules)
	if !f3_n_simp.AmISame(f3_n_corr) {
		TEST.Fatalf("FAIL Div: f3_n_simp != f3_n_corr  ~  %v -> %v  ==  %v", f3_n, f3_n_simp, f3_n_corr)
	} else {
		fmt.Printf("Div:  %v -> %v", f3_n, f3_n_simp)
		TEST.Logf("Div:  %v -> %v", f3_n, f3_n_simp)
	}

	f4_n := parse("(y * x^3) / x^2", vnames)
	f4_n_corr := parse("y * x", vnames)
	f4_n_simp := f4_n.Clone().Simplify(rules)
	if !f4_n_simp.AmISame(f4_n_corr) {
		TEST.Fatalf("FAIL Div: f4_n_simp != f4_n_corr  ~  %v -> %v  ==  %v", f4_n, f4_n_simp, f4_n_corr)
	} else {
		fmt.Printf("Div:  %v -> %v", f4_n, f4_n_simp)
		TEST.Logf("Div:  %v -> %v", f4_n, f4_n_simp)
	}

	f5_n := parse("(y * x^4) / x^2", vnames)
	f5_n_corr := parse("y * x^2", vnames)
	f5_n_simp := f5_n.Clone().Simplify(rules)
	if !f5_n_simp.AmISame(f5_n_corr) {
		TEST.Fatalf("FAIL Div: f5_n_simp != f5_n_corr  ~  %v -> %v  ==  %v", f5_n, f5_n_simp, f5_n_corr)
	} else {
		fmt.Printf("Div:  %v -> %v", f5_n, f5_n_simp)
		TEST.Logf("Div:  %v -> %v", f5_n, f5_n_simp)
	}

	f6_n := parse("(y * x)/ x^2", vnames)
	f6_n_corr := parse("y / x", vnames)
	f6_n_simp := f6_n.Clone().Simplify(rules)
	if !f6_n_simp.AmISame(f6_n_corr) {
		TEST.Fatalf("FAIL Div: f6_n_simp != f6_n_corr  ~  %v -> %v  ==  %v", f6_n, f6_n_simp, f6_n_corr)
	} else {
		fmt.Printf("Div:  %v -> %v", f6_n, f6_n_simp)
		TEST.Logf("Div:  %v -> %v", f6_n, f6_n_simp)
	}

	f7_n := parse("(y * x^2) / x^3", vnames)
	f7_n_corr := parse("y / x", vnames)
	f7_n_simp := f7_n.Clone().Simplify(rules)
	if !f7_n_simp.AmISame(f7_n_corr) {
		TEST.Fatalf("FAIL Div: f7_n_simp != f7_n_corr  ~  %v -> %v  ==  %v", f7_n, f7_n_simp, f7_n_corr)
	} else {
		fmt.Printf("Div:  %v -> %v", f7_n, f7_n_simp)
		TEST.Logf("Div:  %v -> %v", f7_n, f7_n_simp)
	}

	f8_n := parse("(y * x^2) / x^4", vnames)
	f8_n_corr := parse("y / x^2", vnames)
	f8_n_simp := f8_n.Clone().Simplify(rules)
	if !f8_n_simp.AmISame(f8_n_corr) {
		TEST.Fatalf("FAIL Div: f8_n_simp != f8_n_corr  ~  %v -> %v  ==  %v", f8_n, f8_n_simp, f8_n_corr)
	} else {
		fmt.Printf("Div:  %v -> %v", f8_n, f8_n_simp)
		TEST.Logf("Div:  %v -> %v", f8_n, f8_n_simp)
	}

	/*****************************************/
	/*****************************************/

	f1_d := parse("x / (x * y)", vnames)
	f1_d_corr := NewDiv(NewConstant(-1), NewVar(1))
	f1_d_simp := f1_d.Clone().Simplify(rules)
	if !f1_d_simp.AmISame(f1_d_corr) {
		TEST.Fatalf("FAIL Div: f1_d_simp != f1_d_corr  ~  %v -> %v  ==  %v", f1_d, f1_d_simp, f1_d_corr)
	} else {
		fmt.Printf("Div:  %v -> %v", f1_d, f1_d_simp)
		TEST.Logf("Div:  %v -> %v", f1_d, f1_d_simp)
	}

	f2_d := parse("( 2 * x ) / (x * y)", vnames)
	f2_d_corr := parse("2 / y", vnames)
	f2_d_simp := f2_d.Clone().Simplify(rules)
	if !f2_d_simp.AmISame(f2_d_corr) {
		TEST.Fatalf("FAIL Div: f2_d_simp != f2_d_corr  ~  %v -> %v  ==  %v", f2_d, f2_d_simp, f2_d_corr)
	} else {
		fmt.Printf("Div:  %v -> %v", f2_d, f2_d_simp)
		TEST.Logf("Div:  %v -> %v", f2_d, f2_d_simp)
	}

	f3_d := parse("x^2 / (x * y)", vnames)
	f3_d_corr := parse("x / y", vnames)
	f3_d_simp := f3_d.Clone().Simplify(rules)
	if !f3_d_simp.AmISame(f3_d_corr) {
		TEST.Fatalf("FAIL Div: f3_d_simp != f3_d_corr  ~  %v -> %v  ==  %v", f3_d, f3_d_simp, f3_d_corr)
	} else {
		fmt.Printf("Div:  %v -> %v", f3_d, f3_d_simp)
		TEST.Logf("Div:  %v -> %v", f3_d, f3_d_simp)
	}

	f4_d := parse("x^3 / (x^2 * y)", vnames)
	f4_d_corr := parse("x / y", vnames)
	f4_d_simp := f4_d.Clone().Simplify(rules)
	if !f4_d_simp.AmISame(f4_d_corr) {
		TEST.Fatalf("FAIL Div: f4_d_simp != f4_d_corr  ~  %v -> %v  ==  %v", f4_d, f4_d_simp, f4_d_corr)
	} else {
		fmt.Printf("Div:  %v -> %v", f4_d, f4_d_simp)
		TEST.Logf("Div:  %v -> %v", f4_d, f4_d_simp)
	}

	f5_d := parse("x^4 / (x^2 * y)", vnames)
	f5_d_corr := parse("x^2 / y", vnames)
	f5_d_simp := f5_d.Clone().Simplify(rules)
	if !f5_d_simp.AmISame(f5_d_corr) {
		TEST.Fatalf("FAIL Div: f5_d_simp != f5_d_corr  ~  %v -> %v  ==  %v", f5_d, f5_d_simp, f5_d_corr)
	} else {
		fmt.Printf("Div:  %v -> %v", f5_d, f5_d_simp)
		TEST.Logf("Div:  %v -> %v", f5_d, f5_d_simp)
	}

	f6_d := parse("x / (x^2 * y)", vnames)
	f6_d_denom := parse("x * y", vnames)
	f6_d_corr := NewDiv(NewConstant(-1), f6_d_denom)
	f6_d_simp := f6_d.Clone().Simplify(rules)
	if !f6_d_simp.AmISame(f6_d_corr) {
		TEST.Fatalf("FAIL Div: f6_d_simp != f6_d_corr  ~  %v -> %v  ==  %v", f6_d, f6_d_simp, f6_d_corr)
	} else {
		fmt.Printf("Div:  %v -> %v", f6_d, f6_d_simp)
		TEST.Logf("Div:  %v -> %v", f6_d, f6_d_simp)
	}

	f7_d := parse("x^2 / (x^3 * y)", vnames)
	f7_d_denom := parse("x * y", vnames)
	f7_d_corr := NewDiv(NewConstant(-1), f7_d_denom)
	f7_d_simp := f7_d.Clone().Simplify(rules)
	if !f7_d_simp.AmISame(f7_d_corr) {
		TEST.Fatalf("FAIL Div: f7_d_simp != f7_d_corr  ~  %v -> %v  ==  %v", f7_d, f7_d_simp, f7_d_corr)
	} else {
		fmt.Printf("Div:  %v -> %v", f7_d, f7_d_simp)
		TEST.Logf("Div:  %v -> %v", f7_d, f7_d_simp)
	}

	f8_d := parse("x^2 / (x^4 * y)", vnames)
	f8_d_denom := parse("x^2 * y", vnames)
	f8_d_corr := NewDiv(NewConstant(-1), f8_d_denom)
	f8_d_simp := f8_d.Clone().Simplify(rules)
	if !f8_d_simp.AmISame(f8_d_corr) {
		TEST.Fatalf("FAIL Div: f8_d_simp != f8_d_corr  ~  %v -> %v  ==  %v", f8_d, f8_d_simp, f8_d_corr)
	} else {
		fmt.Printf("Div:  %v -> %v", f8_d, f8_d_simp)
		TEST.Logf("Div:  %v -> %v", f8_d, f8_d_simp)
	}

	/*****************************************/
	/*****************************************/

	f1_b := parse("(x * y) / (x * y)", vnames)
	f1_b_corr := NewConstant(-1)
	f1_b_simp := f1_b.Clone().Simplify(rules)
	if !f1_b_simp.AmISame(f1_b_corr) {
		TEST.Fatalf("FAIL Div: f1_b_simp != f1_b_corr  ~  %v -> %v  ==  %v", f1_b, f1_b_simp, f1_b_corr)
	} else {
		fmt.Printf("Div:  %v -> %v", f1_b, f1_b_simp)
		TEST.Logf("Div:  %v -> %v", f1_b, f1_b_simp)
	}

	f2_b := parse("( 2 * x * y ) / (x * y)", vnames)
	f2_b_corr := parse("2", vnames)
	f2_b_simp := f2_b.Clone().Simplify(rules)
	if !f2_b_simp.AmISame(f2_b_corr) {
		TEST.Fatalf("FAIL Div: f2_b_simp != f2_b_corr  ~  %v -> %v  ==  %v", f2_b, f2_b_simp, f2_b_corr)
	} else {
		fmt.Printf("Div:  %v -> %v", f2_b, f2_b_simp)
		TEST.Logf("Div:  %v -> %v", f2_b, f2_b_simp)
	}

	f3_b := parse("(x^2 * y) / (x * y)", vnames)
	f3_b_corr := parse("x", vnames)
	f3_b_simp := f3_b.Clone().Simplify(rules)
	if !f3_b_simp.AmISame(f3_b_corr) {
		TEST.Fatalf("FAIL Div: f3_b_simp != f3_b_corr  ~  %v -> %v  ==  %v", f3_b, f3_b_simp, f3_b_corr)
	} else {
		fmt.Printf("Div:  %v -> %v", f3_b, f3_b_simp)
		TEST.Logf("Div:  %v -> %v", f3_b, f3_b_simp)
	}

	f4_b := parse("(x^3 * y) / (x^2 * y)", vnames)
	f4_b_corr := parse("x", vnames)
	f4_b_simp := f4_b.Clone().Simplify(rules)
	if !f4_b_simp.AmISame(f4_b_corr) {
		TEST.Fatalf("FAIL Div: f4_b_simp != f4_b_corr  ~  %v -> %v  ==  %v", f4_b, f4_b_simp, f4_b_corr)
	} else {
		fmt.Printf("Div:  %v -> %v", f4_b, f4_b_simp)
		TEST.Logf("Div:  %v -> %v", f4_b, f4_b_simp)
	}

	f5_b := parse("(x^4 * y) / (x^2 * y)", vnames)
	f5_b_corr := parse("x^2", vnames)
	f5_b_simp := f5_b.Clone().Simplify(rules)
	if !f5_b_simp.AmISame(f5_b_corr) {
		TEST.Fatalf("FAIL Div: f5_b_simp != f5_b_corr  ~  %v -> %v  ==  %v", f5_b, f5_b_simp, f5_b_corr)
	} else {
		fmt.Printf("Div:  %v -> %v", f5_b, f5_b_simp)
		TEST.Logf("Div:  %v -> %v", f5_b, f5_b_simp)
	}

	f6_b := parse("(x * y) / (x^2 * y)", vnames)
	f6_b_corr := NewDiv(NewConstant(-1), NewVar(0))
	f6_b_simp := f6_b.Clone().Simplify(rules)
	if !f6_b_simp.AmISame(f6_b_corr) {
		TEST.Fatalf("FAIL Div: f6_b_simp != f6_b_corr  ~  %v -> %v  ==  %v", f6_b, f6_b_simp, f6_b_corr)
	} else {
		fmt.Printf("Div:  %v -> %v", f6_b, f6_b_simp)
		TEST.Logf("Div:  %v -> %v", f6_b, f6_b_simp)
	}

	f7_b := parse("(x^2 * y) / (x^3 * y)", vnames)
	f7_b_corr := NewDiv(NewConstant(-1), NewVar(0))
	f7_b_simp := f7_b.Clone().Simplify(rules)
	if !f7_b_simp.AmISame(f7_b_corr) {
		TEST.Fatalf("FAIL Div: f7_b_simp != f7_b_corr  ~  %v -> %v  ==  %v", f7_b, f7_b_simp, f7_b_corr)
	} else {
		fmt.Printf("Div:  %v -> %v", f7_b, f7_b_simp)
		TEST.Logf("Div:  %v -> %v", f7_b, f7_b_simp)
	}

	f8_b := parse("(x^2 * y) / (x^4 * y)", vnames)
	f8_b_corr := NewDiv(NewConstant(-1), NewPowI(NewVar(0), 2))
	f8_b_simp := f8_b.Clone().Simplify(rules)
	if !f8_b_simp.AmISame(f8_b_corr) {
		TEST.Fatalf("FAIL Div: f8_b_simp != f8_b_corr  ~  %v -> %v  ==  %v", f8_b, f8_b_simp, f8_b_corr)
	} else {
		fmt.Printf("Div:  %v -> %v", f8_b, f8_b_simp)
		TEST.Logf("Div:  %v -> %v", f8_b, f8_b_simp)
	}

	/*****************************************/
	/*****************************************/

	f1_B := parse("(x * y) / (x * z)", vnames)
	f1_B_corr := parse("y / z", vnames)
	f1_B_simp := f1_B.Clone().Simplify(rules)
	if !f1_B_simp.AmISame(f1_B_corr) {
		TEST.Fatalf("FAIL Div: f1_B_simp != f1_B_corr  ~  %v -> %v  ==  %v", f1_B, f1_B_simp, f1_B_corr)
	} else {
		fmt.Printf("Div:  %v -> %v", f1_B, f1_B_simp)
		TEST.Logf("Div:  %v -> %v", f1_B, f1_B_simp)
	}

	f2_B := parse("( 2 * x * y) / (x * z)", vnames)
	f2_B_corr := parse("(2 * y) / z", vnames)
	f2_B_simp := f2_B.Clone().Simplify(rules)
	if !f2_B_simp.AmISame(f2_B_corr) {
		TEST.Fatalf("FAIL Div: f2_B_simp != f2_B_corr  ~  %v -> %v  ==  %v", f2_B, f2_B_simp, f2_B_corr)
	} else {
		fmt.Printf("Div:  %v -> %v", f2_B, f2_B_simp)
		TEST.Logf("Div:  %v -> %v", f2_B, f2_B_simp)
	}

	f3_B := parse("(x^2 * y) / (x * z)", vnames)
	f3_B_corr := parse("(x * y) / z", vnames)
	f3_B_simp := f3_B.Clone().Simplify(rules)
	if !f3_B_simp.AmISame(f3_B_corr) {
		TEST.Fatalf("FAIL Div: f3_B_simp != f3_B_corr  ~  %v -> %v  ==  %v", f3_B, f3_B_simp, f3_B_corr)
	} else {
		fmt.Printf("Div:  %v -> %v", f3_B, f3_B_simp)
		TEST.Logf("Div:  %v -> %v", f3_B, f3_B_simp)
	}

	f4_B := parse("(x^3 * y) / (x^2 * z)", vnames)
	f4_B_corr := parse("(x * y) / z", vnames)
	f4_B_simp := f4_B.Clone().Simplify(rules)
	if !f4_B_simp.AmISame(f4_B_corr) {
		TEST.Fatalf("FAIL Div: f4_B_simp != f4_B_corr  ~  %v -> %v  ==  %v", f4_B, f4_B_simp, f4_B_corr)
	} else {
		fmt.Printf("Div:  %v -> %v", f4_B, f4_B_simp)
		TEST.Logf("Div:  %v -> %v", f4_B, f4_B_simp)
	}

	f5_B := parse("(x^4 * y) / (x^2 * z)", vnames)
	f5_B_corr := parse("(x^2 * y) / z", vnames)
	f5_B_simp := f5_B.Clone().Simplify(rules)
	if !f5_B_simp.AmISame(f5_B_corr) {
		TEST.Fatalf("FAIL Div: f5_B_simp != f5_B_corr  ~  %v -> %v  ==  %v", f5_B, f5_B_simp, f5_B_corr)
	} else {
		fmt.Printf("Div:  %v -> %v", f5_B, f5_B_simp)
		TEST.Logf("Div:  %v -> %v", f5_B, f5_B_simp)
	}

	f6_B := parse("(x * y) / (x^2 * z)", vnames)
	f6_B_corr := parse("y / (x * z)", vnames)
	f6_B_simp := f6_B.Clone().Simplify(rules)
	if !f6_B_simp.AmISame(f6_B_corr) {
		TEST.Fatalf("FAIL Div: f6_B_simp != f6_B_corr  ~  %v -> %v  ==  %v", f6_B, f6_B_simp, f6_B_corr)
	} else {
		fmt.Printf("Div:  %v -> %v", f6_B, f6_B_simp)
		TEST.Logf("Div:  %v -> %v", f6_B, f6_B_simp)
	}

	f7_B := parse("(x^2 * y) / (x^3 * z)", vnames)
	f7_B_corr := parse("y / (x * z)", vnames)
	f7_B_simp := f7_B.Clone().Simplify(rules)
	if !f7_B_simp.AmISame(f7_B_corr) {
		TEST.Fatalf("FAIL Div: f7_B_simp != f7_B_corr  ~  %v -> %v  ==  %v", f7_B, f7_B_simp, f7_B_corr)
	} else {
		fmt.Printf("Div:  %v -> %v", f7_B, f7_B_simp)
		TEST.Logf("Div:  %v -> %v", f7_B, f7_B_simp)
	}

	f8_B := parse("(x^2 * y) / (x^4 * z)", vnames)
	f8_B_corr := parse("y / (x^2 * z)", vnames)
	f8_B_simp := f8_B.Clone().Simplify(rules)
	if !f8_B_simp.AmISame(f8_B_corr) {
		TEST.Fatalf("FAIL Div: f8_B_simp != f8_B_corr  ~  %v -> %v  ==  %v", f8_B, f8_B_simp, f8_B_corr)
	} else {
		fmt.Printf("Div:  %v -> %v", f8_B, f8_B_simp)
		TEST.Logf("Div:  %v -> %v", f8_B, f8_B_simp)
	}

	f1_R := parse("( {-1.004803*(x)^2}/{x} + {0.010163*(x)^12}/{x} + 0.995197*x )", vnames)
	f1_R_corr := parse("2 + (3 * x^4)", vnames)
	f1_R_simp := f1_R.Clone().Simplify(rules)
	if !f1_R_simp.AmISame(f1_R_corr) {
		TEST.Fatalf("FAIL Div: f1_R_simp != f1_R_corr  ~  %v -> %v  ==  %v", f1_R, f1_R_simp, f1_R_corr)
	} else {
		fmt.Printf("Div:  %v -> %v", f1_R, f1_R_simp)
		TEST.Logf("Div:  %v -> %v", f1_R, f1_R_simp)
	}

	fmt.Println()
}

// serial := make([]int, 0, 64)
// serial = e.Serial(serial)
// fmt.Printf("\nserial: %v\n", serial)
