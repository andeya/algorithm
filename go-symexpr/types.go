package symexpr

import "sort"

type Leaf struct {
	ExprStats
}

type Unary struct {
	ExprStats
	C Expr
}

type Binary struct {
	ExprStats
	C1, C2 Expr
}

type N_ary struct {
	ExprStats
	CS []Expr
}

// Leaf Nodes
type Time struct {
	Leaf
}

func NewTime() *Time               { return new(Time) }
func (t *Time) ExprType() ExprType { return TIME }
func (t *Time) Clone() Expr        { return NewTime() }

type Var struct {
	Leaf
	P int
}

func NewVar(i int) *Var {
	v := new(Var)
	v.P = i
	return v
}
func (v *Var) ExprType() ExprType { return VAR }
func (v *Var) Clone() Expr        { return NewVar(v.P) }

type Constant struct {
	Leaf
	P int
}

func NewConstant(i int) *Constant {
	c := new(Constant)
	c.P = i
	return c
}
func (c *Constant) ExprType() ExprType { return CONSTANT }
func (c *Constant) Clone() Expr        { return NewConstant(c.P) }

type ConstantF struct {
	Leaf
	F float64
}

func NewConstantF(f float64) *ConstantF {
	c := new(ConstantF)
	c.F = f
	return c
}
func (c *ConstantF) ExprType() ExprType { return CONSTANTF }
func (c *ConstantF) Clone() Expr        { return NewConstantF(c.F) }

type System struct {
	Leaf
	P int
}

func NewSystem(i int) *System {
	s := new(System)
	s.P = i
	return s
}
func (s *System) ExprType() ExprType { return SYSTEM }
func (s *System) Clone() Expr        { return NewSystem(s.P) }

// Unary Operators
type Neg struct {
	Unary
}

func NewNeg(e Expr) *Neg {
	n := new(Neg)
	n.C = e
	return n
}
func (u *Neg) ExprType() ExprType { return NEG }
func (u *Neg) Clone() Expr {
	var C Expr
	if u.C != nil {
		C = u.C.Clone()
	}
	return NewNeg(C)
}

type Abs struct {
	Unary
}

func NewAbs(e Expr) *Abs {
	n := new(Abs)
	n.C = e
	return n
}
func (u *Abs) ExprType() ExprType { return ABS }
func (u *Abs) Clone() Expr {
	var C Expr
	if u.C != nil {
		C = u.C.Clone()
	}
	return NewAbs(C)
}

type Sqrt struct {
	Unary
}

func NewSqrt(e Expr) *Sqrt {
	n := new(Sqrt)
	n.C = e
	return n
}
func (u *Sqrt) ExprType() ExprType { return SQRT }
func (u *Sqrt) Clone() Expr {
	var C Expr
	if u.C != nil {
		C = u.C.Clone()
	}
	return NewSqrt(C)
}

type Sin struct {
	Unary
}

func NewSin(e Expr) *Sin {
	n := new(Sin)
	n.C = e
	return n
}
func (u *Sin) ExprType() ExprType { return SIN }
func (u *Sin) Clone() Expr {
	var C Expr
	if u.C != nil {
		C = u.C.Clone()
	}
	return NewSin(C)
}

type Cos struct {
	Unary
}

func NewCos(e Expr) *Cos {
	n := new(Cos)
	n.C = e
	return n
}
func (u *Cos) ExprType() ExprType { return COS }
func (u *Cos) Clone() Expr {
	var C Expr
	if u.C != nil {
		C = u.C.Clone()
	}
	return NewCos(C)
}

type Tan struct {
	Unary
}

func NewTan(e Expr) *Tan {
	n := new(Tan)
	n.C = e
	return n
}
func (u *Tan) ExprType() ExprType { return TAN }
func (u *Tan) Clone() Expr {
	var C Expr
	if u.C != nil {
		C = u.C.Clone()
	}
	return NewTan(C)
}

type Exp struct {
	Unary
}

func NewExp(e Expr) *Exp {
	n := new(Exp)
	n.C = e
	return n
}
func (u *Exp) ExprType() ExprType { return EXP }
func (u *Exp) Clone() Expr {
	var C Expr
	if u.C != nil {
		C = u.C.Clone()
	}
	return NewExp(C)
}

type Log struct {
	Unary
}

func NewLog(e Expr) *Log {
	n := new(Log)
	n.C = e
	return n
}
func (u *Log) ExprType() ExprType { return LOG }
func (u *Log) Clone() Expr {
	var C Expr
	if u.C != nil {
		C = u.C.Clone()
	}
	return NewLog(C)
}

// Hmmm... Operators
type PowI struct {
	ExprStats
	Base  Expr
	Power int
}

func NewPowI(e Expr, i int) *PowI {
	n := new(PowI)
	if e != nil {
		n.Base = e.Clone()
	}
	n.Power = i
	return n
}
func (u *PowI) ExprType() ExprType { return POWI }
func (u *PowI) Clone() Expr        { return NewPowI(u.Base, u.Power) }

type PowF struct {
	ExprStats
	Base  Expr
	Power float64
}

func NewPowF(b Expr, f float64) *PowF {
	n := new(PowF)
	n.Base = b
	n.Power = f
	return n
}
func (u *PowF) ExprType() ExprType { return POWF }
func (u *PowF) Clone() Expr {
	var base Expr
	if u.Base != nil {
		base = u.Base.Clone()
	}
	return NewPowF(base, u.Power)
}

type PowE struct {
	ExprStats
	Base  Expr
	Power Expr
}

func NewPowE(b, p Expr) *PowE {
	n := new(PowE)
	n.Base = b
	n.Power = p
	return n
}
func (n *PowE) ExprType() ExprType { return POWE }
func (n *PowE) Clone() Expr {
	var base, pow Expr
	if n.Base != nil {
		base = n.Base.Clone()
	}
	if n.Power != nil {
		pow = n.Power.Clone()
	}
	return NewPowE(base, pow)
}

type Div struct {
	ExprStats
	Numer Expr
	Denom Expr
}

func NewDiv(n, d Expr) *Div {
	D := new(Div)
	D.Numer = n
	D.Denom = d
	return D
}
func (n *Div) ExprType() ExprType { return DIV }
func (n *Div) Clone() Expr {
	var N, D Expr
	if n.Numer != nil {
		N = n.Numer.Clone()
	}
	if n.Denom != nil {
		D = n.Denom.Clone()
	}
	return NewDiv(N, D)
}

// N-ary Operators
type Add struct {
	N_ary
}

func NewAdd() *Add {
	a := new(Add)
	a.CS = make([]Expr, 0)
	return a
}
func (n *Add) ExprType() ExprType { return ADD }

func (n *Add) Clone() Expr {
	a := new(Add)
	a.CS = make([]Expr, len(n.CS))
	for i, C := range n.CS {
		if C != nil {
			a.CS[i] = C.Clone()
		}
	}
	return a
}
func (a *Add) Insert(e Expr) {
	if len(a.CS) == cap(a.CS) {
		tmp := make([]Expr, len(a.CS), 2*len(a.CS))
		copy(tmp[:len(a.CS)], a.CS)
		a.CS = tmp
	}
	a.CS = append(a.CS, e)
	sort.Sort(a)
}

type Mul struct {
	N_ary
}

func NewMul() *Mul {
	m := new(Mul)
	m.CS = make([]Expr, 0)
	return m
}
func (n *Mul) ExprType() ExprType { return MUL }
func (n *Mul) Clone() Expr {
	a := NewMul()
	a.CS = make([]Expr, len(n.CS))
	for i, C := range n.CS {
		if C != nil {
			a.CS[i] = C.Clone()
		}
	}
	return a
}
func (a *Mul) Insert(e Expr) {
	if len(a.CS) == cap(a.CS) {
		tmp := make([]Expr, len(a.CS), 2*len(a.CS))
		copy(tmp[:len(a.CS)], a.CS)
		a.CS = tmp
	}
	a.CS = append(a.CS, e)
	sort.Sort(a)
}

func (n *Add) Len() int      { return len(n.CS) }
func (n *Add) Swap(i, j int) { n.CS[i], n.CS[j] = n.CS[j], n.CS[i] }
func (n *Add) Less(i, j int) bool {
	if n.CS[i] == nil && n.CS[j] == nil {
		return false
	}
	if n.CS[i] == nil {
		return false
	}
	if n.CS[j] == nil {
		return true
	}
	return n.CS[i].AmILess(n.CS[j])
}

func (n *Mul) Len() int      { return len(n.CS) }
func (n *Mul) Swap(i, j int) { n.CS[i], n.CS[j] = n.CS[j], n.CS[i] }
func (n *Mul) Less(i, j int) bool {
	if n.CS[i] == nil && n.CS[j] == nil {
		return false
	}
	if n.CS[i] == nil {
		return false
	}
	if n.CS[j] == nil {
		return true
	}
	return n.CS[i].AmILess(n.CS[j])
}
