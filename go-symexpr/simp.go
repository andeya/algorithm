package symexpr

import (
	"math"

	// "fmt"
)

type SimpRules struct {
	GroupAddTerms bool
	ConvertConsts bool
	MulToPow      bool
	MaxPowPow     int

	MulInCoeff int // add must be >= to size, 0 == No

	// state tracking
	InTrig bool
}

func DefaultRules() SimpRules {
	return SimpRules{true, true, true, 12, 0, false}
}

func (l *Leaf) Simplify(rules SimpRules) Expr  { return nil }
func (u *Unary) Simplify(rules SimpRules) Expr { return nil }
func (n *N_ary) Simplify(rules SimpRules) Expr { return nil }

func (n *Time) Simplify(rules SimpRules) Expr { return n }

func (v *Var) Simplify(rules SimpRules) Expr { return v }

func (c *Constant) Simplify(rules SimpRules) Expr { return c }

func (c *ConstantF) Simplify(rules SimpRules) Expr {
	if math.IsNaN(c.F) || math.IsInf(c.F, 0) {
		return nil
	}
	// "close" to system value ??? TODO

	return c
}

func (s *System) Simplify(rules SimpRules) Expr { return s }

func (u *Neg) Simplify(rules SimpRules) (ret Expr) {
	if u.C == nil {
		return nil
	}
	u.C = u.C.Simplify(rules)
	if u.C == nil {
		return nil
	}
	switch u.C.ExprType() {
	case NULL:
		ret = nil
		u.C = nil
	case NEG:
		ret = u.C
		u.C = nil
	case CONSTANT:
		ret = u.C
		u.C = nil
	case CONSTANTF:
		c := u.C.(*ConstantF)
		c.F *= -1.0
		ret = u.C
		u.C = nil
	case MUL:
		m := u.C.(*Mul)
		if m.CS[0].ExprType() == CONSTANT { // Constants should always be first operand in Mul
			ret = u.C
			u.C = nil
		} else if m.CS[0].ExprType() == CONSTANTF { // Constants should always be first operand in Mul
			c := m.CS[0].(*ConstantF)
			c.F *= -1.0
			ret = u.C
			u.C = nil
		}
	default: // no simplification
		ret = u
	}
	return ret
}

func (u *Abs) Simplify(rules SimpRules) (ret Expr) {
	if u.C == nil {
		return nil
	}
	u.C = u.C.Simplify(rules)
	if u.C == nil {
		return nil
	}
	switch u.C.ExprType() {
	case NULL:
		ret = nil
		u.C = nil
	case ABS:
		ret = u.C
		u.C = nil
	default: // no simplification
		ret = u
	}
	return ret
}

func (u *Sqrt) Simplify(rules SimpRules) (ret Expr) {
	if u.C == nil {
		return nil
	}
	u.C = u.C.Simplify(rules)
	if u.C == nil {
		return nil
	}
	switch u.C.ExprType() {
	case NULL:
		ret = nil
		u.C = nil
	case SQRT:
		s := u.C.(*Sqrt)
		ret = NewPowF(s.C, 0.25)
		u.C = nil
	case CONSTANT:
		ret = u.C
		u.C = nil
	case CONSTANTF:
		ret = u.C
		u.C = nil
		ret.(*ConstantF).F = math.Sqrt(ret.(*ConstantF).F)
	case POWF:
		p := u.C.(*PowF)
		p.Power *= 0.5
		ret = p
		u.C = nil
	case POWI:
		p := u.C.(*PowI)
		if p.Power == 2 {
			ret = p.Base
			p.Base = nil
		} else if p.Power%2 == 0 {
			p.Power /= 2
			ret = p
			u.C = nil
		} else {
			f := NewPowF(p.Base, float64(p.Power)*0.5)
			ret = f
			p.Base = nil
		}

	default: // no simplification
		ret = u
	}
	return ret
}

func (u *Sin) Simplify(rules SimpRules) (ret Expr) {
	// fmt.Printf("Simp Sin:  %v\n", u.C)
	if !rules.InTrig {
		rules.InTrig = true
	} else {
		tmp := u.C
		u.C = nil
		return tmp
	}
	if u.C == nil {
		return nil
	}
	u.C = u.C.Simplify(rules)
	if u.C == nil {
		return nil
	}
	// fmt.Printf("      %v\n", u.C)
	switch u.C.ExprType() {
	case NULL, SIN, COS, TAN, CONSTANT:
		ret = u.C
		u.C = nil
	case CONSTANTF:
		ret = u.C
		u.C = nil
		ret.(*ConstantF).F = math.Sin(ret.(*ConstantF).F)
	default: // no simplification
		ret = u
	}
	return ret
}

func (u *Cos) Simplify(rules SimpRules) (ret Expr) {
	if !rules.InTrig {
		rules.InTrig = true
	} else {
		tmp := u.C
		u.C = nil
		return tmp
	}
	if u.C == nil {
		return nil
	}
	u.C = u.C.Simplify(rules)
	if u.C == nil {
		return nil
	}
	switch u.C.ExprType() {
	case NULL, SIN, COS, TAN, CONSTANT:
		ret = u.C
		u.C = nil
	case CONSTANTF:
		ret = u.C
		u.C = nil
		ret.(*ConstantF).F = math.Cos(ret.(*ConstantF).F)
	default: // no simplification
		ret = u
	}
	return ret
}

func (u *Tan) Simplify(rules SimpRules) (ret Expr) {
	if !rules.InTrig {
		rules.InTrig = true
	} else {
		tmp := u.C
		u.C = nil
		return tmp
	}
	if u.C == nil {
		return nil
	}
	u.C = u.C.Simplify(rules)
	if u.C == nil {
		return nil
	}
	switch u.C.ExprType() {
	case NULL, SIN, COS, TAN, CONSTANT:
		ret = u.C
		u.C = nil
	case CONSTANTF:
		ret = u.C
		u.C = nil
		ret.(*ConstantF).F = math.Tan(ret.(*ConstantF).F)
	default: // no simplification
		ret = u
	}
	return ret
}

func (u *Exp) Simplify(rules SimpRules) (ret Expr) {
	if u.C == nil {
		return nil
	}
	u.C = u.C.Simplify(rules)
	if u.C == nil {
		return nil
	}
	switch u.C.ExprType() {
	case NULL, CONSTANT:
		ret = u.C
		u.C = nil
	case CONSTANTF:
		ret = u.C
		u.C = nil
		ret.(*ConstantF).F = math.Exp(ret.(*ConstantF).F)
	default: // no simplification
		ret = u
	}
	return ret
}

func (u *Log) Simplify(rules SimpRules) (ret Expr) {
	if u.C == nil {
		return nil
	}
	u.C = u.C.Simplify(rules)
	if u.C == nil {
		return nil
	}
	switch u.C.ExprType() {
	case NULL, CONSTANT:
		ret = u.C
		u.C = nil
	case CONSTANTF:
		ret = u.C
		u.C = nil
		ret.(*ConstantF).F = math.Log(ret.(*ConstantF).F)
	default: // no simplification
		ret = u
	}
	return ret
}

func (u *PowI) Simplify(rules SimpRules) Expr {
	var (
		ret Expr     = u
		t   ExprType = NULL
	)
	if u.Base != nil {
		// serial := make([]int, 0, 64)
		// serial = u.Base.Serial(serial)
		// fmt.Printf("PowI-presort:  %v   %v\n", u.Base, serial)

		u.Base = u.Base.Simplify(rules)
		if u.Base == nil {
			return nil
		}

		// serial2 := make([]int, 0, 64)
		// serial2 = u.Base.Serial(serial2)
		// fmt.Printf("PowI-postsort:  %v   %v\n", u.Base, serial2)
		t = u.Base.ExprType()
	}

	if u.Power > rules.MaxPowPow {
		u.Power = rules.MaxPowPow
	} else if u.Power < -rules.MaxPowPow {
		u.Power = -rules.MaxPowPow
	}

	if u.Power == 0 {
		ret = &ConstantF{F: 1}
	} else if u.Power == 1 {
		ret = u.Base
		u.Base = nil
	} else {
		switch t {
		case NULL, CONSTANT:
			ret = u.Base
			u.Base = nil
		case CONSTANTF:
			ret = u.Base
			u.Base = nil
			ret.(*ConstantF).F = math.Pow(ret.(*ConstantF).F, float64(u.Power))
		case MUL:
			if u.Base.HasConst() {
				m := u.Base.(*Mul)
				mret := NewMul()
				mret.Insert(m.CS[0])
				m.CS[0] = nil
				u.Base = m.Simplify(rules)
				tmp := NewPowF(u.Base, float64(u.Power))
				mret.Insert(tmp)
				ret = mret
			}
		case POWI:
			p := u.Base.(*PowI)
			u.Base = p.Base
			u.Power = u.Power * p.Power
			if u.Power > rules.MaxPowPow {
				u.Power = rules.MaxPowPow
			} else if u.Power < -rules.MaxPowPow {
				u.Power = -rules.MaxPowPow
			}
			// ret = NewPowF(u.Base, float64(u.Power))
			ret = u
		default: // no simplification
			// ret = NewPowF(u.Base, float64(u.Power))
			ret = u
		}
	}
	// serial := make([]int, 0, 64)
	// serial = ret.Serial(serial)
	// fmt.Printf("PowI-ret:  %v  %v\n", ret, serial)
	return ret
}

func (u *PowF) Simplify(rules SimpRules) Expr {
	var (
		ret Expr     = u
		t   ExprType = NULL
	)
	if u.Base != nil {
		u.Base = u.Base.Simplify(rules)
		if u.Base == nil {
			return nil
		}
		t = u.Base.ExprType()
	}

	if u.Power == 0 {
		ret = &ConstantF{F: 1}
	} else if u.Power == 1 {
		ret = u.Base
		u.Base = nil
	} else {
		switch t {
		case NULL, CONSTANT:
			ret = u.Base
			u.Base = nil
		case CONSTANTF:
			ret = u.Base
			u.Base = nil
			ret.(*ConstantF).F = math.Pow(ret.(*ConstantF).F, float64(u.Power))
		default: // no simplification
			ret = NewPowI(u.Base, int(u.Power))
		}
	}
	return ret
}

func (n *PowE) Simplify(rules SimpRules) Expr {
	var (
		ret    Expr     = n
		t1, t2 ExprType = NULL, NULL
	)
	if n.Base != nil {
		n.Base = n.Base.Simplify(rules)
		t1 = n.Base.ExprType()
	}
	if n.Power != nil {
		n.Power = n.Power.Simplify(rules)
		t2 = n.Power.ExprType()
	}
	if t1 == NULL && t2 == NULL {
		return &Null{}
	} else if t1 == NULL {
		ret = n.Power
		n.Base = nil
		n.Power = nil
	} else if t2 == NULL {
		ret = n.Base
		n.Base = nil
		n.Power = nil
	} else if n.Base.ExprType() == n.Power.ExprType() &&
		n.Base.ExprType() == CONSTANTF {
		ret = n.Base
		ret.(*ConstantF).F = math.Pow(ret.(*ConstantF).F, n.Power.(*ConstantF).F)
		n.Base = nil
		n.Power = nil
	}
	return ret
}

func (n *Div) Simplify(rules SimpRules) Expr {
	var (
		ret    Expr     = n
		t1, t2 ExprType = NULL, NULL
	)
	if n.Numer != nil {
		n.Numer = n.Numer.Simplify(rules)
		if n.Numer != nil {
			t1 = n.Numer.ExprType()
		}
	}
	if n.Denom != nil {
		n.Denom = n.Denom.Simplify(rules)
		if n.Denom != nil {
			t2 = n.Denom.ExprType()
		}
	}

	// check to see if nulls/nils
	if t1 == NULL && t2 == NULL {
		ret = nil
	} else if t1 == NULL {
		ret = n.Denom
		n.Numer = nil
		n.Denom = nil
		return ret
	} else if t2 == NULL {
		ret = n.Numer
		n.Numer = nil
		n.Denom = nil
		return ret
	}

	// check to see if divs
	if t1 == DIV && t2 == DIV {
		//  a/b // c/d => ad/bc
		ndiv := n.Numer.(*Div)
		ddiv := n.Denom.(*Div)
		nmul := NewMul()
		nmul.Insert(ndiv.Numer)
		nmul.Insert(ddiv.Denom)
		dmul := NewMul()
		dmul.Insert(ndiv.Denom)
		dmul.Insert(ddiv.Numer)
		ndiv.Numer = nil
		ndiv.Denom = nil
		ddiv.Numer = nil
		ddiv.Denom = nil
		n.Numer = nmul.Simplify(rules)
		n.Denom = dmul.Simplify(rules)
		t1 = n.Numer.ExprType()
		t2 = n.Denom.ExprType()

	} else if t1 == DIV {
		//  a/b // c   => a/bc
		div := n.Numer.(*Div)
		n.Numer = div.Numer
		mul := NewMul()
		mul.Insert(div.Denom)
		mul.Insert(n.Denom)
		div.Numer = nil
		div.Denom = nil
		n.Denom = mul.Simplify(rules)
		t1 = n.Numer.ExprType()
		t2 = n.Denom.ExprType()

	} else if t2 == DIV {
		//  a // c/d   => ad/c
		div := n.Denom.(*Div)
		n.Denom = div.Numer
		mul := NewMul()
		mul.Insert(n.Numer)
		mul.Insert(div.Denom)
		div.Numer = nil
		div.Denom = nil
		n.Numer = mul.Simplify(rules)
		t1 = n.Numer.ExprType()
		t2 = n.Denom.ExprType()
	}

	// cancel / simp like terms
	if t1 == t2 && t1 == CONSTANT {
		ret = n.Numer
		n.Numer = nil
		n.Denom = nil
		return ret
	} else if t1 == CONSTANT && t2 == MUL {
		d := n.Denom.(*Mul)
		if d.CS[0].ExprType() == CONSTANT {
			d.CS[0] = nil
			n.Denom = d.Simplify(rules)
		}
		return ret
	} else if t2 == CONSTANT && t1 == MUL {
		d := n.Numer.(*Mul)
		if d.CS[0].ExprType() != CONSTANT {
			d.Insert(NewConstant(-1))
		}
		ret = n.Numer
		n.Numer = nil
		n.Denom = nil
		return ret
	} else if t1 == t2 && t1 == CONSTANTF {
		c := n.Numer.(*ConstantF)
		c.F /= n.Denom.(*ConstantF).F
		ret = n.Numer
		n.Numer = nil
		n.Denom = nil
		return ret
	} else if t1 == CONSTANTF && t2 == MUL {
		d := n.Denom.(*Mul)
		if d.CS[0].ExprType() == CONSTANTF {
			c := n.Numer.(*ConstantF)
			c.F /= d.CS[0].(*ConstantF).F
			d.CS[0] = nil
			n.Denom = d.Simplify(rules)
		}
		return ret
	} else if t2 == CONSTANTF && t1 == MUL {
		d := n.Numer.(*Mul)
		if d.CS[0].ExprType() != CONSTANTF {
			c := NewConstantF(1.0)
			c.F /= n.Denom.(*ConstantF).F
			d.Insert(c)
		} else {
			c := d.CS[0].(*ConstantF)
			c.F /= n.Denom.(*ConstantF).F
		}
		ret = n.Numer
		n.Numer = nil
		n.Denom = nil
		return ret
	} else {
		ret = n.cancelLikeTerms(rules)
	}

	return ret
}

func (D *Div) cancelLikeTerms(rules SimpRules) Expr {

	t1 := D.Numer.ExprType()
	t2 := D.Denom.ExprType()
	var ret Expr = D

	// fmt.Println("cancel", t1, t2)

	if D.Numer.AmISame(D.Denom) {
		return NewConstantF(1.0)
	}

	changed := false
	if t1 == t2 && t1 == MUL {
		mn := D.Numer.(*Mul)
		md := D.Denom.(*Mul)
		for i, n := range mn.CS {
			if n == nil {
				continue
			}
			for j, d := range md.CS {
				if d == nil {
					continue
				}
				// fmt.Println("cancel", i, n.ExprType(), j, d.ExprType())

				if n.AmISame(d) {
					mn.CS[i] = nil
					md.CS[j] = nil
					changed = true
					break
				}
				if n.ExprType() == POWI && d.ExprType() == POWI {
					np, dp := n.(*PowI), d.(*PowI)
					if np.Base.AmISame(dp.Base) {
						if np.Power > dp.Power {
							np.Power -= dp.Power
							md.CS[j] = nil
							changed = true
							continue
						} else if np.Power < dp.Power {
							dp.Power -= np.Power
							mn.CS[i] = nil
							changed = true
							break
						} else { // same and should cancel
							mn.CS[i] = nil
							md.CS[j] = nil
							changed = true
							break
						}
					}
				} else if n.ExprType() == POWI {
					np := n.(*PowI)
					if np.Base.AmISame(d) {
						np.Power -= 1
						md.CS[j] = nil
						changed = true
						continue
					}
				} else if d.ExprType() == POWI {
					nd := d.(*PowI)
					if nd.Base.AmISame(n) {
						nd.Power -= 1
						mn.CS[i] = nil
						changed = true
						break
					}
				}

			}

		}
	} else if t1 == MUL {
		if D.Numer.AmIAlmostSame(D.Denom) {
			return NewConstantF(1.0)
		}
		mn := D.Numer.(*Mul)
		d := D.Denom
		for i, n := range mn.CS {
			if n.AmISame(d) {
				mn.CS[i] = nil
				D.Denom = nil
				changed = true
				break
			}
			if n.ExprType() == POWI && d.ExprType() == POWI {
				np, dp := n.(*PowI), d.(*PowI)
				if np.Base.AmISame(dp.Base) {
					if np.Power > dp.Power {
						np.Power -= dp.Power
						D.Denom = nil
						changed = true
						break
					} else if np.Power < dp.Power {
						dp.Power -= np.Power
						mn.CS[i] = nil
						changed = true
						break
					} else { // same and should cancel
						// will we even get here?
						mn.CS[i] = nil
						D.Denom = nil
						changed = true
						break
					}
				}
			} else if n.ExprType() == POWI {
				np := n.(*PowI)
				if np.Base.AmISame(d) {
					np.Power -= 1
					D.Denom = nil
					changed = true
					break
				}
			} else if d.ExprType() == POWI {
				nd := d.(*PowI)
				if nd.Base.AmISame(n) {
					nd.Power -= 1
					mn.CS[i] = nil
					changed = true
					break
				}
			}
		}
	} else if t2 == MUL {
		if D.Denom.AmIAlmostSame(D.Numer) {
			return NewConstantF(1.0)
		}
		n := D.Numer
		md := D.Denom.(*Mul)
		for i, d := range md.CS {
			if n.AmISame(d) {
				D.Numer = NewConstantF(1.0)
				md.CS[i] = nil
				changed = true
				break
			}
			if n.ExprType() == POWI && d.ExprType() == POWI {
				np, dp := n.(*PowI), d.(*PowI)
				if np.Base.AmISame(dp.Base) {
					if np.Power > dp.Power {
						np.Power -= dp.Power
						md.CS[i] = nil
						changed = true
						break
					} else if np.Power < dp.Power {
						dp.Power -= np.Power
						D.Numer = NewConstantF(1.0)
						changed = true
						break
					} else { // same and should cancel
						// will we even get here?
						D.Numer = nil
						md.CS[i] = nil
						changed = true
						break
					}
				}
			} else if n.ExprType() == POWI {
				np := n.(*PowI)
				if np.Base.AmISame(d) {
					np.Power -= 1
					md.CS[i] = nil
					changed = true
					break
				}
			} else if d.ExprType() == POWI {
				nd := d.(*PowI)
				if nd.Base.AmISame(n) {
					nd.Power -= 1
					D.Numer = NewConstantF(1.0)
					changed = true
					break
				}
			}
		}
	} else {
		n, d := D.Numer, D.Denom
		if n.AmISame(d) {
			D.Numer = NewConstantF(1.0)
			D.Denom = nil
			changed = true
		} else if n.ExprType() == POWI && d.ExprType() == POWI {
			// fmt.Println("Got Here 1")
			np, dp := n.(*PowI), d.(*PowI)
			if np.Base.AmISame(dp.Base) {
				if np.Power > dp.Power {
					np.Power -= dp.Power
					D.Denom = nil
					changed = true
				} else if np.Power < dp.Power {
					dp.Power -= np.Power
					D.Numer = NewConstantF(1.0)
					changed = true
				} else { // same and should cancel
					// will we even get here?
					D.Numer = NewConstantF(1.0)
					D.Denom = nil
					changed = true
				}
			}
		} else if n.ExprType() == POWI {
			// fmt.Println("Got Here 2")
			np := n.(*PowI)
			if np.Base.AmISame(d) {
				np.Power -= 1
				D.Denom = nil
				changed = true
			}
		} else if d.ExprType() == POWI {
			// fmt.Println("Got Here 3")
			nd := d.(*PowI)
			if nd.Base.AmISame(n) {
				nd.Power -= 1
				D.Numer = NewConstantF(1.0)
				changed = true
			}
		}
	}

	if changed {
		t1, t2 = NULL, NULL
		if D.Numer != nil {
			D.Numer = D.Numer.Simplify(rules)
			if D.Numer != nil {
				t1 = D.Numer.ExprType()
			}
		}
		if D.Denom != nil {
			D.Denom = D.Denom.Simplify(rules)
			if D.Denom != nil {
				t2 = D.Denom.ExprType()
			}
		}

		// check to see if nulls/nils
		if t1 == NULL && t2 == NULL {
			ret = nil
		} else if t1 == NULL {
			D.Numer = NewConstantF(1.0)
			ret = D
			return ret
		} else if t2 == NULL {
			ret = D.Numer
			D.Numer = nil
			D.Denom = nil
			return ret
		}
	}

	return ret
}

func (n *Add) Simplify(rules SimpRules) Expr {
	var ret Expr = n
	hasSome := false
	// fmt.Printf("n.CS = %v\n",n.CS)
	for i, C := range n.CS {
		if C != nil {
			n.CS[i] = C.Simplify(rules)
			if n.CS[i] == nil {
				continue
			}
			if n.CS[i].ExprType() == NULL {
				n.CS[i] = nil
			} else {
				hasSome = true
			}
		}
	}

	if !hasSome {
		return nil
	}

	changed := true
	for changed {
		// fmt.Printf("iter = %v\n",n.CS)
		// fmt.Printf("ADD:  %v\n", n)
		changed = false
		if rules.ConvertConsts {
			changed = changed || gatherAddTerms(n)
			// fmt.Printf("GTH:  %v\n", n)
			changed = changed || groupAddTerms(n)
			// fmt.Printf("GRP:  %v\n", n)
		} else {
			changed = changed || gatherAddTerms(n)
			// fmt.Printf("GTH:  %v\n", n)
			changed = changed || groupAddTermsF(n)
			// fmt.Printf("GRP:  %v\n", n)
		}
		n.Sort()
		// fmt.Printf("SRT:  %v\n", n)
	}
	// fmt.Printf("POST:  %v\n\n\n", n)

	cnt := countTerms(n.CS)
	if cnt == 0 {
		ret = nil
	}
	if cnt == 1 {
		ret = n.CS[0]
		n.CS[0] = nil
	}
	// fmt.Printf("exit = %v\n",n.CS)

	return ret
}

func (n *Mul) Simplify(rules SimpRules) Expr {

	var ret Expr = n

	changed := true
	for changed {
		for i, C := range n.CS {
			if C != nil {
				n.CS[i] = C.Simplify(rules)
				if n.CS[i] == nil {
					continue
				}
				if n.CS[i].ExprType() == NULL {
					n.CS[i] = nil
				}
			}
		}

		// fmt.Printf("MUL:  %v\n", n)
		changed = false
		if rules.ConvertConsts {
			changed = changed || gatherMulTerms(n)
			// fmt.Printf("GTH:  %v\n", n)
			changed = changed || groupMulTerms(n)
			// fmt.Printf("GRP:  %v\n", n)
		} else {
			changed = changed || gatherMulTermsF(n)
			// fmt.Printf("GTH:  %v\n", n)
			changed = changed || groupMulTermsF(n)
			// fmt.Printf("GRP:  %v\n", n)
		}
		n.Sort()
		// fmt.Printf("SRT:  %v\n", n)
	}
	// fmt.Printf("POST:  %v\n\n\n", n)

	cnt := countTerms(n.CS)
	if cnt == 0 {
		ret = nil
	}
	if cnt == 1 {
		// fmt.Printf("ARRRRGGGGG   %v\n", n.CS)
		ret = n.CS[0]
		n.CS[0] = nil
	}
	return ret
}

func countTerms(terms []Expr) int {
	cnt := 0
	for _, e := range terms {
		if e != nil {
			cnt++
		}
	}
	return cnt
}

// this function left aligns children terms ( ie move nils to end of terms[] )
// and returns the number of children
func leftAlignTerms(terms []Expr) int {
	cnt, nilp := 0, -1
	for i, e := range terms {
		if e != nil {
			// fmt.Printf( "TERMS(%d/%d): %v\n", i,nilp, terms )
			cnt++
			if nilp >= 0 {
				terms[nilp], terms[i] = terms[i], nil
				nilp++
				// find next nil spot
				for nilp < len(terms) {
					if terms[nilp] == nil {
						break
					} else {
						nilp++
					}
				}
				if nilp >= len(terms) {
					break
				} // the outer loop
			}
		} else if nilp < 0 {
			nilp = i
		}
	}
	return cnt
}

// pulls recursive adds into the current add
func gatherAddTerms(n *Add) (changed bool) {
	terms := make([]Expr, 0)
	// terms := n.CS
	// fmt.Printf("n.CS: %v\n", n.CS[:])
	for i, e := range n.CS {
		if e == nil {
			continue
		}
		if e.ExprType() == ADD {
			// fmt.Printf("got here\n")
			changed = true
			a := e.(*Add)
			for j, E := range a.CS {
				if E == nil {
					continue
				}
				terms = append(terms, E)
				a.CS[j] = nil
			}
			// rem := leftAlignTerms(a.CS[:])
			// if rem == 0 {
			n.CS[i] = nil
			// }
			// leftAlignTerms(terms)
		} else {
			terms = append(terms, e)
		}
	}
	// fmt.Printf("terms: %v\n", terms)
	n.CS = terms
	return changed
}

// pulls recursive muls into the current mul
func gatherMulTerms(n *Mul) (changed bool) {
	terms := make([]Expr, 0)
	hasCoeff := false
	var coeff *ConstantF = nil
	hasDiv := -1
	for i, e := range n.CS {
		if e == nil {
			continue
		}
		switch e.ExprType() {

		case CONSTANT:
			if !hasCoeff {
				terms = append(terms, e)
				hasCoeff = true
			} else {
				changed = true
			}

		case CONSTANTF:
			// changed = true
			if coeff == nil {
				// hasCoeff = true
				coeff = e.(*ConstantF)
				terms = append(terms, e)
				// terms = append(terms, NewConstant(-1))
			} else {
				coeff.F *= e.(*ConstantF).F
			}
			n.CS[i] = nil

		case NEG:
			changed = true
			neg := e.(*Neg)
			if !hasCoeff {
				mul := NewMul()
				mul.Insert(NewConstantF(-1))
				mul.Insert(neg.C)
				neg.C = nil
				terms = append(terms, mul)
				hasCoeff = true
			} else {
				terms = append(terms, neg.C)
				neg.C = nil
			}

		case DIV:
			if hasDiv == -1 {
				hasDiv = i // store the first div occurance
			}
			terms = append(terms, e)
		case MUL:
			changed = true
			a := e.(*Mul)
			leftAlignTerms(a.CS[:])
			for j, E := range a.CS {
				if E == nil {
					continue
				}
				if e.ExprType() == CONSTANT {
					if !hasCoeff {
						hasCoeff = true
						terms = append(terms, E)
					}
					a.CS[i] = nil
				} else {
					terms = append(terms, E)
					a.CS[j] = nil
				}
			}
			rem := leftAlignTerms(a.CS[:])
			if rem == 0 {
				n.CS[i] = nil
			}
			leftAlignTerms(terms)

		default:
			terms = append(terms, e)
		}
	}
	n.CS = terms

	if hasDiv > -1 && len(n.CS) > 1 {
		changed = true
		numer := NewMul()
		denom := NewMul()
		for i, e := range n.CS {
			if e == nil {
				continue
			}
			switch e.ExprType() {
			case DIV:
				d := e.(*Div)
				numer.Insert(d.Numer)
				denom.Insert(d.Denom)
			// case POWI:
			// case POWF:
			default:
				numer.Insert(e)
			}
			n.CS[i] = nil
		}
		n.Insert(NewDiv(numer, denom))
	}
	return changed
}

// group like terms in an addition
func groupAddTerms(n *Add) (changed bool) {
	terms := n.CS
	for i, x := range terms {
		// fmt.Printf("TERMS:  %v\n", terms)
		if x == nil {
			continue
		}
		xT := x.ExprType()
		same := false
		for j := i + 1; j < len(terms); j++ {
			y := terms[j]
			// fmt.Printf("%d-%d  %v  %v\n", i, j, x, y)
			if y == nil {
				continue
			}

			//  F(x) +  F(x)
			// cF(x) + dF(x)
			// cF(x) +  F(x)
			//  F(x) + cF(x)
			if x.AmIAlmostSame(y) || y.AmIAlmostSame(x) { // both directions ensures Mul.AmIAlmostSame()... with only one FAIL on F(x) + cF(x)
				same = true
				terms[j] = nil
			}

			//  F(x) -  F(x)
			// cF(x) - dF(x)
			// cF(x) -  F(x)
			//  F(x) - cF(x)
			if xT == NEG {
				c := x.(*Neg).C
				if y.AmIAlmostSame(c) {
					same = true
					terms[j] = nil
				}
			}
			yT := y.ExprType()
			if yT == NEG {
				c := y.(*Neg).C
				if x.AmIAlmostSame(c) {
					same = true
					terms[j] = nil
				}
			}

		}
		if same {
			changed = true
			// fmt.Printf("groupAddTerms: %v\n", terms)
			// extract x.C if neg
			if xT == NEG {
				x = x.(*Neg).C
				xT = x.ExprType()
			}
			if xT == MUL {
				mul := x.(*Mul)
				if !(mul.CS[0].ExprType() == CONSTANT) {
					mul.Insert(NewConstantF(1.0))
				}
			} else {
				mul := NewMul()
				mul.Insert(NewConstantF(1.0))
				mul.Insert(x)
				terms[i] = mul
			}
		}

	}
	return
}

// group like terms in a multiplication
func groupMulTerms(m *Mul) (changed bool) {
	terms := m.CS
	L := len(terms)
	hasConst := false
	var coeff *ConstantF = nil

	for i, x := range terms {
		if x == nil {
			continue
		}

		// setting up current term
		powSum := 1.0

		xT := x.ExprType()
		if xT == CONSTANT && !hasConst {
			hasConst = true
		}

		if xT == CONSTANTF {
			if coeff == nil {
				coeff = x.(*ConstantF)
			} else {
				coeff.F *= x.(*ConstantF).F
				terms[i] = nil
				continue
			}
		}

		var xC Expr
		xTb := NULL
		switch xT {
		case NEG:
			xC = x.(*Neg).C
			xTb = xC.ExprType()
		case POWI:
			p := x.(*PowI)
			xC = p.Base
			xTb = p.Base.ExprType()
			powSum = float64(p.Power)
		case POWF:
			p := x.(*PowF)
			xC = p.Base
			xTb = p.Base.ExprType()
			powSum = p.Power
		case POWE:
			p := x.(*PowE)
			if p.Power.ExprType() == CONSTANTF {
				xC = p.Base
				xTb = p.Base.ExprType()
				powSum = p.Power.(*ConstantF).F
			}
		}

		for j := i + 1; j < L; j++ {
			y := terms[j]
			if y == nil {
				continue
			}

			// setting up following term(s)
			powY := 1.0
			yT := y.ExprType()
			if yT == CONSTANT && !hasConst {
				hasConst = true
			}

			if yT == CONSTANTF {
				if coeff == nil {
					coeff = y.(*ConstantF)
				} else {
					coeff.F *= y.(*ConstantF).F
					terms[j] = nil
					continue
				}
			}

			var yC Expr
			yTb := NULL
			switch yT {
			case NEG:
				yC = y.(*Neg).C
				yTb = yC.ExprType()
			case POWI:
				p := y.(*PowI)
				yC = p.Base
				yTb = p.Base.ExprType()
				powY = float64(p.Power)
			case POWF:
				p := y.(*PowF)
				yC = p.Base
				yTb = p.Base.ExprType()
				powY = p.Power
			case POWE:
				p := y.(*PowE)
				if p.Power.ExprType() == CONSTANTF {
					yC = p.Base
					yTb = p.Base.ExprType()
					powY = p.Power.(*ConstantF).F
				}
			}

			// fmt.Println("Types: ", xT, xTb, yT, yTb)

			// This is the actual comparison Code
			same := false
			if xT == yT {
				if x.AmISame(y) {
					same = true
					powSum += powY
				}
				if xTb == yTb && xTb != NULL {
					if xC.AmISame(yC) {
						same = true
						powSum += powY
					}
				}
			} else if xTb == yT {
				if xC.AmISame(y) {
					same = true
					powSum += powY
				}
			} else if xT == yTb {
				if x.AmISame(yC) {
					same = true
					powSum += powY
				}
			} else if xTb == yTb && xTb != NULL {
				if xC.AmISame(yC) {
					same = true
					powSum += powY
				}
			}

			// Check the results of camparison update the terms
			if same {
				// fmt.Printf("GotHERE\n")
				terms[j] = nil
				changed = true
				if powSum == 0 {
					terms[i] = nil // remove the lhs term
				} else if powSum == 1 {
					// what is lhs?
					// if PowI || Neg || ? :: then extract Base
					if xTb != NULL {
						switch xT {
						case NEG:
							terms[i] = x.(*Neg).C
						case POWI:
							terms[i] = x.(*PowI).Base
						case POWF:
							terms[i] = x.(*PowF).Base
						}
					}
				} else {
					// whole or fractional power?
					flr := math.Floor(powSum)
					dim := math.Dim(powSum, flr)
					if dim == 0 {
						// whole power
						base := x
						if xT != POWI {
							if xTb != NULL {
								switch xT {
								case NEG:
									base = x.(*Neg).C
								case POWF:
									base = x.(*PowF).Base
								}
							}
							base = NewPowI(base, int(powSum))
						} else {
							base.(*PowI).Power = int(powSum)
						}
						terms[i] = base
					} else {
						// fractional power
						base := x
						if xT != POWF {
							if xTb != NULL {
								switch xT {
								case NEG:
									base = x.(*Neg).C
								case POWI:
									base = x.(*PowI).Base
								}
							}
							base = NewPowF(base, powSum)
						} else {
							base.(*PowF).Power = powSum
						}
					}
				}
			}
		}

	}
	return
}

// pulls recursive muls into the current mul
func gatherMulTermsF(n *Mul) (changed bool) {
	terms := make([]Expr, 0)
	hasCoeff := false
	var coeff *ConstantF = nil
	hasDiv := -1
	for i, e := range n.CS {
		if e == nil {
			continue
		}
		switch e.ExprType() {

		case CONSTANT:
			if !hasCoeff {
				terms = append(terms, e)
				hasCoeff = true
			} else {
				changed = true
			}

		case CONSTANTF:
			// changed = true
			if coeff == nil {
				// hasCoeff = true
				coeff = e.(*ConstantF)
				terms = append(terms, e)
				// terms = append(terms, NewConstant(-1))
			} else {
				coeff.F *= e.(*ConstantF).F
			}
			n.CS[i] = nil

		case NEG:
			changed = true
			neg := e.(*Neg)
			if !hasCoeff {
				mul := NewMul()
				mul.Insert(NewConstantF(-1))
				mul.Insert(neg.C)
				neg.C = nil
				terms = append(terms, mul)
				hasCoeff = true
			} else {
				terms = append(terms, neg.C)
				neg.C = nil
			}

		case DIV:
			if hasDiv == -1 {
				hasDiv = i // store the first div occurance
			}
			terms = append(terms, e)
		case MUL:
			changed = true
			a := e.(*Mul)
			leftAlignTerms(a.CS[:])
			for j, E := range a.CS {
				if E == nil {
					continue
				}
				if e.ExprType() == CONSTANT || e.ExprType() == CONSTANTF {
					if !hasCoeff {
						hasCoeff = true
						terms = append(terms, E)
					}
					a.CS[i] = nil
				} else {
					terms = append(terms, E)
					a.CS[j] = nil
				}
			}
			rem := leftAlignTerms(a.CS[:])
			if rem == 0 {
				n.CS[i] = nil
			}
			leftAlignTerms(terms)

		default:
			terms = append(terms, e)
		}
	}
	n.CS = terms

	if hasDiv > -1 && len(n.CS) > 1 {
		changed = true
		numer := NewMul()
		denom := NewMul()
		for i, e := range n.CS {
			if e == nil {
				continue
			}
			switch e.ExprType() {
			case DIV:
				d := e.(*Div)
				numer.Insert(d.Numer)
				denom.Insert(d.Denom)
			// case POWI:
			// case POWF:
			default:
				numer.Insert(e)
			}
			n.CS[i] = nil
		}
		n.Insert(NewDiv(numer, denom))
	}
	return changed
}

// group like terms in an addition
func groupAddTermsF(n *Add) (changed bool) {
	terms := n.CS
	for i, x := range terms {
		// fmt.Printf("TERMS:  %v\n", terms)
		if x == nil {
			continue
		}
		xT := x.ExprType()
		same := false
		sum := 1.0
		if xT == MUL && x.(*Mul).CS[0].ExprType() == CONSTANTF {
			sum = x.(*Mul).CS[0].(*ConstantF).F
		} else if xT == NEG {
			sum = -1.0
		}
		for j := i + 1; j < len(terms); j++ {
			y := terms[j]
			// fmt.Printf("%d-%d  %v  %v\n", i, j, x, y)
			if y == nil {
				continue
			}
			yT := y.ExprType()

			//  F(x) +  F(x)
			// cF(x) + dF(x)
			// cF(x) +  F(x)
			//  F(x) + cF(x)
			if x.AmIAlmostSame(y) || y.AmIAlmostSame(x) { // both directions ensures Mul.AmIAlmostSame()... with only one FAIL on F(x) + cF(x)
				if yT == MUL && y.(*Mul).CS[0].ExprType() == CONSTANTF {
					sum += y.(*Mul).CS[0].(*ConstantF).F
				} else {
					sum += 1.0
				}
				same = true
				terms[j] = nil
			}

			//  F(x) -  F(x)
			// cF(x) - dF(x)
			// cF(x) -  F(x)
			//  F(x) - cF(x)
			if xT == NEG {
				// fmt.Printf("GOT HERE xT NEG\n")
				c := x.(*Neg).C
				if y.AmIAlmostSame(c) {
					// fmt.Printf("GOT HERE xT NEG y almostSame\n")
					if yT == MUL && y.(*Mul).CS[0].ExprType() == CONSTANTF {
						sum += y.(*Mul).CS[0].(*ConstantF).F
					} else {
						sum += 1.0
					}
					same = true
					terms[j] = nil
				}
			}
			if yT == NEG {
				c := y.(*Neg).C
				if x.AmIAlmostSame(c) {
					if c.ExprType() == MUL && c.(*Mul).CS[0].ExprType() == CONSTANTF {
						sum -= c.(*Mul).CS[0].(*ConstantF).F
					} else {
						sum -= 1.0
					}
					same = true
					terms[j] = nil
				}
			}

		}
		if same {
			changed = true
			// fmt.Printf("sum = %v\n", sum)
			// fmt.Printf("groupAddTerms: %v\n", terms)
			if sum == 0.0 {
				terms[i] = NewConstantF(sum)
				continue
			} else if sum == 1.0 {
				if xT == NEG {
					terms[i] = x.(*Neg).C
				}
				continue
			} else if sum == -1.0 {
				if xT != NEG {
					terms[i] = NewNeg(x)
				}
				continue
			} else {

				// extract x.C if neg
				if xT == NEG {
					// fmt.Printf("xT is NEG\n")
					x = x.(*Neg).C
					xT = x.ExprType()
				}
				if xT == MUL {
					// fmt.Printf("xT is MUL\n")
					mul := x.(*Mul)
					if !(mul.CS[0].ExprType() == CONSTANTF) {
						mul.Insert(NewConstantF(sum))
					} else {
						mul.CS[0].(*ConstantF).F = sum
					}
				} else {
					mul := NewMul()
					mul.Insert(NewConstantF(sum))
					mul.Insert(x)
					terms[i] = mul
				}
			}
		}

	}
	return
}

// group like terms in a multiplication
func groupMulTermsF(m *Mul) (changed bool) {
	terms := m.CS
	L := len(terms)
	hasConst := false
	var coeff *ConstantF = nil

	for i, x := range terms {
		if x == nil {
			continue
		}

		// setting up current term
		powSum := 1.0

		xT := x.ExprType()
		if xT == CONSTANT && !hasConst {
			hasConst = true
		}

		if xT == CONSTANTF {
			if coeff == nil {
				coeff = x.(*ConstantF)
			} else {
				coeff.F *= x.(*ConstantF).F
				terms[i] = nil
				continue
			}
		}

		var xC Expr
		xTb := NULL
		switch xT {
		case NEG:
			xC = x.(*Neg).C
			xTb = xC.ExprType()
		case POWI:
			p := x.(*PowI)
			xC = p.Base
			xTb = p.Base.ExprType()
			powSum = float64(p.Power)
		case POWF:
			p := x.(*PowF)
			xC = p.Base
			xTb = p.Base.ExprType()
			powSum = p.Power
		}

		for j := i + 1; j < L; j++ {
			y := terms[j]
			if y == nil {
				continue
			}

			// setting up following term(s)
			powY := 1.0
			yT := y.ExprType()
			if yT == CONSTANT && !hasConst {
				hasConst = true
			}

			if yT == CONSTANTF {
				if coeff == nil {
					coeff = y.(*ConstantF)
				} else {
					coeff.F *= y.(*ConstantF).F
					terms[j] = nil
					continue
				}
			}

			var yC Expr
			yTb := NULL
			switch yT {
			case NEG:
				yC = y.(*Neg).C
				yTb = yC.ExprType()
			case POWI:
				p := y.(*PowI)
				yC = p.Base
				yTb = p.Base.ExprType()
				powY = float64(p.Power)
			case POWF:
				p := y.(*PowF)
				yC = p.Base
				yTb = p.Base.ExprType()
				powY = p.Power
			}

			// This is the actual comparison Code
			same := false
			if xT == yT {
				if x.AmISame(y) {
					same = true
					powSum += powY
				}
			} else if xTb == yT {
				if xC.AmISame(y) {
					same = true
					powSum += powY
				}
			} else if xT == yTb {
				if x.AmISame(yC) {
					same = true
					powSum += powY
				}
			} else if xTb == yTb && xTb != NULL {
				if xC.AmISame(yC) {
					same = true
					powSum += powY
				}
			}

			// Check the results of camparison update the terms
			if same {
				terms[j] = nil
				changed = true
				if powSum == 0 {
					terms[i] = nil // remove the lhs term
				} else if powSum == 1 {
					// what is lhs?
					// if PowI || Neg || ? :: then extract Base
					if xTb != NULL {
						switch xT {
						case NEG:
							terms[i] = x.(*Neg).C
						case POWI:
							terms[i] = x.(*PowI).Base
						case POWF:
							terms[i] = x.(*PowF).Base
						}
					}
				} else {
					// // whole or fractional power?
					// flr := math.Floor(powSum)
					// dim := math.Dim(powSum, flr)
					// if dim == 0 {
					// 	// whole power
					// 	base := x
					// 	if xT != POWI {
					// 		if xTb != NULL {
					// 			switch xT {
					// 			case NEG:
					// 				base = x.(*Neg).C
					// 			case POWF:
					// 				base = x.(*PowF).Base
					// 			}
					// 		}
					// 		base = NewPowI(base, int(powSum))
					// 	} else {
					// 		base.(*PowI).Power = int(powSum)
					// 	}
					// 	terms[i] = base
					// } else {
					// fractional power
					base := x
					if xT != POWF {
						if xTb != NULL {
							switch xT {
							case NEG:
								base = x.(*Neg).C
							case POWI:
								base = x.(*PowI).Base
							}
						}
						base = NewPowF(base, powSum)
					} else {
						base.(*PowF).Power = powSum
					}
					// }
				}
			}
		}

	}
	return
}
