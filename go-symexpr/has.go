package symexpr

func (l *Leaf) HasVar() bool         { return false }
func (l *Leaf) HasVarI(i int) bool   { return false }
func (l *Leaf) NumVar() int          { return 0 }
func (l *Leaf) HasConst() bool       { return false }
func (l *Leaf) HasConstI(i int) bool { return false }
func (l *Leaf) NumConstants() int    { return 0 }

func (u *Unary) HasVar() bool         { return u.C.HasVar() }
func (u *Unary) HasVarI(i int) bool   { return u.C.HasVarI(i) }
func (u *Unary) NumVar() int          { return u.C.NumVar() }
func (u *Unary) HasConst() bool       { return u.C.HasConst() }
func (u *Unary) HasConstI(i int) bool { return u.C.HasConstI(i) }
func (u *Unary) NumConstants() int    { return u.C.NumConstants() }

func (n *N_ary) HasVar() bool {
	for _, C := range n.CS {
		if C != nil && C.HasVar() {
			return true
		}
	}
	return false
}
func (n *N_ary) HasVarI(i int) bool {
	for _, C := range n.CS {
		if C != nil && C.HasVarI(i) {
			return true
		}
	}
	return false
}
func (n *N_ary) NumVar() int {
	sum := 0
	for _, C := range n.CS {
		if C != nil {
			sum += C.NumVar()
		}
	}
	return sum
}
func (n *N_ary) HasConst() bool {
	for _, C := range n.CS {
		if C != nil && C.HasConst() {
			return true
		}
	}
	return false
}
func (n *N_ary) HasConstI(i int) bool {
	for _, C := range n.CS {
		if C != nil && C.HasConstI(i) {
			return true
		}
	}
	return false
}
func (n *N_ary) NumConstants() int {
	sum := 0
	for _, C := range n.CS {
		if C != nil {
			sum += C.NumConstants()
		}
	}
	return sum
}

// Overload specifics

func (v *Var) HasVar() bool       { return true }
func (v *Var) HasVarI(i int) bool { return v.P == i }
func (v *Var) NumVar() int        { return 1 }

func (c *Constant) HasConst() bool       { return true }
func (c *Constant) HasConstI(i int) bool { return c.P == i }
func (c *Constant) NumConstants() int    { return 1 }

func (c *ConstantF) HasConst() bool       { return false }
func (c *ConstantF) HasConstI(i int) bool { return false }

// func (c *ConstantF) NumConstants() int      { return 1 }   NOT COUNTING THESE

func (u *PowI) HasVar() bool         { return u.Base.HasVar() }
func (u *PowI) HasVarI(i int) bool   { return u.Base.HasVarI(i) }
func (u *PowI) NumVar() int          { return u.Base.NumVar() }
func (u *PowI) HasConst() bool       { return u.Base.HasConst() }
func (u *PowI) HasConstI(i int) bool { return u.Base.HasConstI(i) }
func (u *PowI) NumConstants() int    { return u.Base.NumConstants() }

func (u *PowF) HasVar() bool         { return u.Base.HasVar() }
func (u *PowF) HasVarI(i int) bool   { return u.Base.HasVarI(i) }
func (u *PowF) NumVar() int          { return u.Base.NumVar() }
func (u *PowF) HasConst() bool       { return u.Base.HasConst() }
func (u *PowF) HasConstI(i int) bool { return u.Base.HasConstI(i) }
func (u *PowF) NumConstants() int    { return u.Base.NumConstants() }

func (n *PowE) HasVar() bool         { return n.Base.HasVar() || n.Power.HasVar() }
func (n *PowE) HasVarI(i int) bool   { return n.Base.HasVarI(i) || n.Power.HasVarI(i) }
func (n *PowE) NumVar() int          { return n.Base.NumVar() + n.Power.NumVar() }
func (n *PowE) HasConst() bool       { return n.Base.HasConst() || n.Power.HasConst() }
func (n *PowE) HasConstI(i int) bool { return n.Base.HasConstI(i) || n.Power.HasConstI(i) }
func (n *PowE) NumConstants() int    { return n.Base.NumConstants() + n.Power.NumConstants() }

func (n *Div) HasVar() bool {
	ret := false
	if n.Numer != nil {
		ret = ret || n.Numer.HasVar()
	}
	if n.Denom != nil {
		ret = ret || n.Denom.HasVar()
	}
	return ret
}
func (n *Div) HasVarI(i int) bool {
	ret := false
	if n.Numer != nil {
		ret = ret || n.Numer.HasVarI(i)
	}
	if n.Denom != nil {
		ret = ret || n.Denom.HasVarI(i)
	}
	return ret
}
func (n *Div) NumVar() int {
	ret := 0
	if n.Numer != nil {
		ret = ret + n.Numer.NumVar()
	}
	if n.Denom != nil {
		ret = ret + n.Denom.NumVar()
	}
	return ret
}
func (n *Div) HasConst() bool {
	ret := false
	if n.Numer != nil {
		ret = ret || n.Numer.HasConst()
	}
	if n.Denom != nil {
		ret = ret || n.Denom.HasConst()
	}
	return ret
}
func (n *Div) HasConstI(i int) bool {
	ret := false
	if n.Numer != nil {
		ret = ret || n.Numer.HasConstI(i)
	}
	if n.Denom != nil {
		ret = ret || n.Denom.HasConstI(i)
	}
	return ret
}
func (n *Div) NumConstants() int {
	ret := 0
	if n.Numer != nil {
		ret = ret + n.Numer.NumConstants()
	}
	if n.Denom != nil {
		ret = ret + n.Denom.NumConstants()
	}
	return ret
}
