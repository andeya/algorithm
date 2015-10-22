package symexpr

import (
	"sort"

	// "fmt"
)

func (n *Time) AmILess(r Expr) bool       { return TIME < r.ExprType() }
func (n *Time) AmIEqual(r Expr) bool      { return r.ExprType() == TIME }
func (n *Time) AmISame(r Expr) bool       { return r.ExprType() == TIME }
func (n *Time) AmIAlmostSame(r Expr) bool { return r.ExprType() == TIME }
func (n *Time) Sort()                     { return }

func (v *Var) AmILess(r Expr) bool {
	if VAR < r.ExprType() {
		return true
	}
	if VAR > r.ExprType() {
		return false
	}
	return v.P < r.(*Var).P
}
func (v *Var) AmIEqual(r Expr) bool      { return r.ExprType() == VAR && r.(*Var).P == v.P }
func (v *Var) AmISame(r Expr) bool       { return r.ExprType() == VAR && r.(*Var).P == v.P }
func (v *Var) AmIAlmostSame(r Expr) bool { return r.ExprType() == VAR && r.(*Var).P == v.P }
func (v *Var) Sort()                     { return }

func (c *Constant) AmILess(r Expr) bool {
	if CONSTANT < r.ExprType() {
		return true
	}
	if CONSTANT > r.ExprType() {
		return false
	}
	return c.P < r.(*Constant).P
}
func (c *Constant) AmIEqual(r Expr) bool      { return r.ExprType() == CONSTANT && r.(*Constant).P == c.P }
func (c *Constant) AmISame(r Expr) bool       { return r.ExprType() == CONSTANT }
func (c *Constant) AmIAlmostSame(r Expr) bool { return r.ExprType() == CONSTANT }
func (c *Constant) Sort()                     { return }

func (c *ConstantF) AmILess(r Expr) bool {
	if CONSTANTF < r.ExprType() {
		return true
	}
	if CONSTANTF > r.ExprType() {
		return false
	}
	return c.F < r.(*ConstantF).F
}
func (c *ConstantF) AmIEqual(r Expr) bool {
	return r.ExprType() == CONSTANTF && r.(*ConstantF).F == c.F
}
func (c *ConstantF) AmISame(r Expr) bool       { return r.ExprType() == CONSTANTF }
func (c *ConstantF) AmIAlmostSame(r Expr) bool { return r.ExprType() == CONSTANTF }
func (c *ConstantF) Sort()                     { return }

func (s *System) AmILess(r Expr) bool {
	if SYSTEM < r.ExprType() {
		return true
	}
	if SYSTEM > r.ExprType() {
		return false
	}
	return s.P < r.(*System).P
}
func (s *System) AmIEqual(r Expr) bool      { return r.ExprType() == SYSTEM && r.(*System).P == s.P }
func (s *System) AmISame(r Expr) bool       { return r.ExprType() == SYSTEM && r.(*System).P == s.P }
func (s *System) AmIAlmostSame(r Expr) bool { return r.ExprType() == SYSTEM && r.(*System).P == s.P }
func (s *System) Sort()                     { return }

func (u *Neg) AmILess(r Expr) bool {
	if NEG < r.ExprType() {
		return true
	}
	if NEG > r.ExprType() {
		return false
	}
	return u.C.AmILess(r.(*Neg).C)
}
func (u *Neg) AmIEqual(r Expr) bool { return r.ExprType() == NEG && u.C.AmIEqual(r.(*Neg).C) }
func (u *Neg) AmISame(r Expr) bool  { return r.ExprType() == NEG && u.C.AmISame(r.(*Neg).C) }
func (u *Neg) AmIAlmostSame(r Expr) bool {
	return r.ExprType() == NEG && u.C.AmIAlmostSame(r.(*Neg).C)
}
func (u *Neg) Sort() { u.C.Sort() }

func (u *Abs) AmILess(r Expr) bool {
	if ABS < r.ExprType() {
		return true
	}
	if ABS > r.ExprType() {
		return false
	}
	return u.C.AmILess(r.(*Abs).C)
}
func (u *Abs) AmIEqual(r Expr) bool { return r.ExprType() == ABS && u.C.AmIEqual(r.(*Abs).C) }
func (u *Abs) AmISame(r Expr) bool  { return r.ExprType() == ABS && u.C.AmISame(r.(*Abs).C) }
func (u *Abs) AmIAlmostSame(r Expr) bool {
	return r.ExprType() == ABS && u.C.AmIAlmostSame(r.(*Abs).C)
}
func (u *Abs) Sort() { u.C.Sort() }

func (u *Sqrt) AmILess(r Expr) bool {
	if SQRT < r.ExprType() {
		return true
	}
	if SQRT > r.ExprType() {
		return false
	}
	return u.C.AmILess(r.(*Sqrt).C)
}
func (u *Sqrt) AmIEqual(r Expr) bool { return r.ExprType() == SQRT && u.C.AmIEqual(r.(*Sqrt).C) }
func (u *Sqrt) AmISame(r Expr) bool  { return r.ExprType() == SQRT && u.C.AmISame(r.(*Sqrt).C) }
func (u *Sqrt) AmIAlmostSame(r Expr) bool {
	return r.ExprType() == SQRT && u.C.AmIAlmostSame(r.(*Sqrt).C)
}
func (u *Sqrt) Sort() { u.C.Sort() }

func (u *Sin) AmILess(r Expr) bool {
	if SIN < r.ExprType() {
		return true
	}
	if SIN > r.ExprType() {
		return false
	}
	return u.C.AmILess(r.(*Sin).C)
}
func (u *Sin) AmIEqual(r Expr) bool { return r.ExprType() == SIN && u.C.AmIEqual(r.(*Sin).C) }
func (u *Sin) AmISame(r Expr) bool  { return r.ExprType() == SIN && u.C.AmISame(r.(*Sin).C) }
func (u *Sin) AmIAlmostSame(r Expr) bool {
	return r.ExprType() == SIN && u.C.AmIAlmostSame(r.(*Sin).C)
}
func (u *Sin) Sort() { u.C.Sort() }

func (u *Cos) AmILess(r Expr) bool {
	if COS < r.ExprType() {
		return true
	}
	if COS > r.ExprType() {
		return false
	}
	return u.C.AmILess(r.(*Cos).C)
}
func (u *Cos) AmIEqual(r Expr) bool { return r.ExprType() == COS && u.C.AmIEqual(r.(*Cos).C) }
func (u *Cos) AmISame(r Expr) bool  { return r.ExprType() == COS && u.C.AmISame(r.(*Cos).C) }
func (u *Cos) AmIAlmostSame(r Expr) bool {
	return r.ExprType() == COS && u.C.AmIAlmostSame(r.(*Cos).C)
}
func (u *Cos) Sort() { u.C.Sort() }

func (u *Tan) AmILess(r Expr) bool {
	if TAN < r.ExprType() {
		return true
	}
	if TAN > r.ExprType() {
		return false
	}
	return u.C.AmILess(r.(*Tan).C)
}
func (u *Tan) AmIEqual(r Expr) bool { return r.ExprType() == TAN && u.C.AmIEqual(r.(*Tan).C) }
func (u *Tan) AmISame(r Expr) bool  { return r.ExprType() == TAN && u.C.AmISame(r.(*Tan).C) }
func (u *Tan) AmIAlmostSame(r Expr) bool {
	return r.ExprType() == TAN && u.C.AmIAlmostSame(r.(*Tan).C)
}
func (u *Tan) Sort() { u.C.Sort() }

func (u *Exp) AmILess(r Expr) bool {
	if EXP < r.ExprType() {
		return true
	}
	if EXP > r.ExprType() {
		return false
	}
	return u.C.AmILess(r.(*Exp).C)
}
func (u *Exp) AmIEqual(r Expr) bool { return r.ExprType() == EXP && u.C.AmIEqual(r.(*Exp).C) }
func (u *Exp) AmISame(r Expr) bool  { return r.ExprType() == EXP && u.C.AmISame(r.(*Exp).C) }
func (u *Exp) AmIAlmostSame(r Expr) bool {
	return r.ExprType() == EXP && u.C.AmIAlmostSame(r.(*Exp).C)
}
func (u *Exp) Sort() { u.C.Sort() }

func (u *Log) AmILess(r Expr) bool {
	if LOG < r.ExprType() {
		return true
	}
	if LOG > r.ExprType() {
		return false
	}
	return u.C.AmILess(r.(*Log).C)
}
func (u *Log) AmIEqual(r Expr) bool { return r.ExprType() == LOG && u.C.AmIEqual(r.(*Log).C) }
func (u *Log) AmISame(r Expr) bool  { return r.ExprType() == LOG && u.C.AmISame(r.(*Log).C) }
func (u *Log) AmIAlmostSame(r Expr) bool {
	return r.ExprType() == LOG && u.C.AmIAlmostSame(r.(*Log).C)
}
func (u *Log) Sort() { u.C.Sort() }

func (u *PowI) AmILess(r Expr) bool {
	if POWI < r.ExprType() {
		return true
	}
	if POWI > r.ExprType() {
		return false
	}
	if u.Base.AmILess(r.(*PowI).Base) {
		return true
	}
	if r.(*PowI).Base.AmILess(u.Base) {
		return false
	}
	return u.Power < r.(*PowI).Power
}
func (u *PowI) AmIEqual(r Expr) bool {
	return r.ExprType() == POWI && r.(*PowI).Power == u.Power && u.Base.AmIEqual(r.(*PowI).Base)
}
func (u *PowI) AmISame(r Expr) bool {
	return r.ExprType() == POWI && r.(*PowI).Power == u.Power && u.Base.AmISame(r.(*PowI).Base)
}
func (u *PowI) AmIAlmostSame(r Expr) bool {
	return r.ExprType() == POWI && r.(*PowI).Power == u.Power && u.Base.AmIAlmostSame(r.(*PowI).Base)
}
func (u *PowI) Sort() { u.Base.Sort() }

func (u *PowF) AmILess(r Expr) bool {
	if POWF < r.ExprType() {
		return true
	}
	if POWF > r.ExprType() {
		return false
	}
	if u.Base.AmILess(r.(*PowF).Base) {
		return true
	}
	if r.(*PowF).Base.AmILess(u.Base) {
		return false
	}
	return u.Power < r.(*PowF).Power
}
func (u *PowF) AmIEqual(r Expr) bool {
	return r.ExprType() == POWF && r.(*PowF).Power == u.Power && u.Base.AmIEqual(r.(*PowF).Base)
}
func (u *PowF) AmISame(r Expr) bool {
	return r.ExprType() == POWF && r.(*PowF).Power == u.Power && u.Base.AmISame(r.(*PowF).Base)
}
func (u *PowF) AmIAlmostSame(r Expr) bool {
	return r.ExprType() == POWF && r.(*PowF).Power == u.Power && u.Base.AmIAlmostSame(r.(*PowF).Base)
}
func (u *PowF) Sort() { u.Base.Sort() }

func (n *PowE) AmILess(r Expr) bool {
	if POWE < r.ExprType() {
		return true
	}
	if POWE > r.ExprType() {
		return false
	}
	if n.Base.AmILess(r.(*PowE).Base) {
		return true
	}
	if r.(*PowE).Base.AmILess(n.Base) {
		return false
	}
	return n.Power.AmILess(r.(*PowE).Power)
}
func (n *PowE) AmIEqual(r Expr) bool {
	return r.ExprType() == POWE && n.Base.AmIEqual(r.(*PowE).Base) && n.Power.AmIEqual(r.(*PowE).Power)
}
func (n *PowE) AmISame(r Expr) bool {
	return r.ExprType() == POWE && n.Base.AmISame(r.(*PowE).Base) && n.Power.AmISame(r.(*PowE).Power)
}
func (n *PowE) AmIAlmostSame(r Expr) bool {
	return r.ExprType() == POWE && n.Base.AmIAlmostSame(r.(*PowE).Base) && n.Power.AmIAlmostSame(r.(*PowE).Power)
}
func (n *PowE) Sort() { n.Base.Sort(); n.Power.Sort() }

func (n *Div) AmILess(r Expr) bool {
	if r == nil {
		return false
	}
	if DIV < r.ExprType() {
		return true
	}
	if DIV > r.ExprType() {
		return false
	}
	rp := r.(*Div)
	if n.Numer.AmILess(rp.Numer) {
		return true
	}
	if rp.Numer.AmILess(n.Numer) {
		return false
	}
	return n.Denom.AmILess(r.(*Div).Denom)
}
func (n *Div) AmIEqual(r Expr) bool {
	if r == nil || r.ExprType() != DIV {
		return false
	}
	rp := r.(*Div)
	if (n.Numer != nil && rp.Numer == nil) || (n.Numer == nil && rp.Numer != nil) {
		return false
	}
	if (n.Denom != nil && rp.Denom == nil) || (n.Denom == nil && rp.Denom != nil) {
		return false
	}
	return r.ExprType() == DIV && n.Numer.AmIEqual(rp.Numer) && rp.Denom.AmIEqual(rp.Denom)
}
func (n *Div) AmISame(r Expr) bool {
	if r == nil || r.ExprType() != DIV {
		return false
	}
	rp := r.(*Div)
	if (n.Numer != nil && rp.Numer == nil) || (n.Numer == nil && rp.Numer != nil) {
		return false
	}
	if (n.Denom != nil && rp.Denom == nil) || (n.Denom == nil && rp.Denom != nil) {
		return false
	}
	return r.ExprType() == DIV && n.Numer.AmISame(rp.Numer) && n.Denom.AmISame(rp.Denom)
}
func (n *Div) AmIAlmostSame(r Expr) bool {
	if r == nil || r.ExprType() != DIV {
		return false
	}
	rp := r.(*Div)
	if (n.Numer != nil && rp.Numer == nil) || (n.Numer == nil && rp.Numer != nil) {
		return false
	}
	if (n.Denom != nil && rp.Denom == nil) || (n.Denom == nil && rp.Denom != nil) {
		return false
	}
	return r.ExprType() == DIV && n.Numer.AmIAlmostSame(rp.Numer) && n.Denom.AmIAlmostSame(rp.Denom)
}
func (n *Div) Sort() {
	if n.Numer != nil {
		n.Numer.Sort()
	}
	if n.Denom != nil {
		n.Denom.Sort()
	}
}

func (n *Add) AmILess(r Expr) bool {
	if ADD < r.ExprType() {
		return true
	}
	if ADD > r.ExprType() {
		return false
	}
	m := r.(*Add)
	ln := len(n.CS)
	lm := len(m.CS)
	if ln < lm {
		return true
	}
	if lm < ln {
		return false
	}
	for i, C := range n.CS {
		if C.AmILess(m.CS[i]) {
			return true
		}
		if m.CS[i].AmILess(C) {
			return false
		}
	}
	return false
}
func (n *Add) AmIEqual(r Expr) bool {
	if r.ExprType() != ADD {
		return false
	}
	m := r.(*Add)
	if len(n.CS) != len(m.CS) {
		return false
	}
	for i, C := range n.CS {
		if !C.AmIEqual(m.CS[i]) {
			return false
		}
		//     if m.CS[i].AmILess( C ) { return false }
	}
	return true
}
func (n *Add) AmISame(r Expr) bool {
	if r.ExprType() != ADD {
		return false
	}
	m := r.(*Add)
	if len(n.CS) != len(m.CS) {
		return false
	}
	for i, C := range n.CS {
		if !C.AmISame(m.CS[i]) {
			return false
		}
		//     if m.CS[i].AmILess( C ) { return false }
	}
	return true
}
func (n *Add) AmIAlmostSame(r Expr) bool {
	if r.ExprType() != ADD {
		return false
	}
	m := r.(*Add)
	if len(n.CS) != len(m.CS) {
		return false
	}
	same := true
	for i, C := range n.CS {
		if !C.AmIAlmostSame(m.CS[i]) {
			return false
		}
	}
	return same
}
func (n *Add) Sort() {
	for _, C := range n.CS {
		if C != nil {
			C.Sort()
		}
	}
	sort.Sort(ExprArray(n.CS))
	i := len(n.CS) - 1
	for i >= 0 && n.CS[i] == nil {
		n.CS = n.CS[:i]
		i = len(n.CS) - 1
	}
	return
}

func (n *Mul) AmILess(r Expr) bool {
	if r.ExprType() != MUL {
		if len(n.CS) == 2 {
			if n.CS[0].ExprType() == CONSTANT && n.CS[1].AmISame(r) {
				return true
			}
			if n.CS[0].ExprType() == CONSTANTF && n.CS[1].AmISame(r) {
				return true
			}
		}
		return false
	}
	m := r.(*Mul)
	ln, lm := len(n.CS), len(m.CS)
	sn, sm := 0, 0
	if n.CS[0].ExprType() == CONSTANT || n.CS[0].ExprType() == CONSTANTF {
		sn++
	}
	if m.CS[0].ExprType() == CONSTANT || m.CS[0].ExprType() == CONSTANTF {
		sm++
	}

	if ln-sn != lm-sm {
		return ln-sn < lm-sm
	}
	for i, j := sn, sm; i < ln && j < lm; {
		if n.CS[i].AmILess(m.CS[j]) {
			return true
		}
		if m.CS[j].AmILess(n.CS[i]) {
			return false
		}
		i++
		j++
	}
	return false
}

func (n *Mul) AmIEqual(r Expr) bool {
	if r.ExprType() != MUL {
		return false
	}
	m := r.(*Mul)
	if len(n.CS) != len(m.CS) {
		return false
	}
	for i, C := range n.CS {
		if !C.AmIEqual(m.CS[i]) {
			return false
		}
		//     if m.CS[i].AmILess( C ) { return false }
	}
	return true
}

func (n *Mul) AmISame(r Expr) bool {
	if r.ExprType() != MUL {
		return false
	}
	m := r.(*Mul)
	if len(n.CS) != len(m.CS) {
		return false
	}
	for i, C := range n.CS {
		if !C.AmISame(m.CS[i]) {
			return false
		}
		//     if m.CS[i].AmILess( C ) { return false }
	}
	return true
}

func (n *Mul) AmIAlmostSame(r Expr) bool {
	if r.ExprType() != MUL {
		// fmt.Printf("~MUL:  %v  %v\n", n, r)
		if len(n.CS) == 2 {
			if n.CS[0].ExprType() == CONSTANT && n.CS[1].AmIAlmostSame(r) {
				return true
			}
			if n.CS[0].ExprType() == CONSTANTF && n.CS[1].AmIAlmostSame(r) {
				return true
			}
		}
		return false
	}
	// fmt.Printf("MUL:  %v  %v\n", n, r)

	m := r.(*Mul)
	ln, lm := len(n.CS), len(m.CS)
	sn, sm := 0, 0
	if n.CS[0].ExprType() == CONSTANT || n.CS[0].ExprType() == CONSTANTF {
		sn++
	}
	if m.CS[0].ExprType() == CONSTANT || m.CS[0].ExprType() == CONSTANTF {
		sm++
	}
	// fmt.Printf("lens: %d %d %d %d\n", sn, ln, sm, lm)

	if ln-sn != lm-sm {
		return false
	}
	for i, j := sn, sm; i < ln && j < lm; {
		// fmt.Printf("COMPARE:  %v   %v", n.CS[i], m.CS[j])
		if !n.CS[i].AmIAlmostSame(m.CS[j]) {
			// fmt.Println("SAME")
			return false
		}
		i++
		j++
	}
	// fmt.Printf("\n")
	return true
}

func (n *Mul) Sort() {
	for _, C := range n.CS {
		if C != nil {
			C.Sort()
		}
	}
	sort.Sort(ExprArray(n.CS))
	i := len(n.CS) - 1
	for i >= 0 && n.CS[i] == nil {
		n.CS = n.CS[:i]
		i = len(n.CS) - 1
	}
	return
}
