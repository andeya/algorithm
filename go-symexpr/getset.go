package symexpr

func (n *Time) GetExpr(pos *int) Expr {
	if (*pos) == 0 {
		return n
	}
	(*pos)--
	return nil
}
func (n *Time) SetExpr(pos *int, e Expr) (replace_me, replaced bool) {
	if (*pos) == 0 {
		return true, false
	}
	(*pos)--
	return false, false
}

func (v *Var) GetExpr(pos *int) Expr {
	if (*pos) == 0 {
		return v
	}
	(*pos)--
	return nil
}
func (v *Var) SetExpr(pos *int, e Expr) (replace_me, replaced bool) {
	if (*pos) == 0 {
		return true, false
	}
	(*pos)--
	return false, false
}

func (c *Constant) GetExpr(pos *int) Expr {
	if (*pos) == 0 {
		return c
	}
	(*pos)--
	return nil
}
func (c *Constant) SetExpr(pos *int, e Expr) (replace_me, replaced bool) {
	if (*pos) == 0 {
		return true, false
	}
	(*pos)--
	return false, false
}

func (c *ConstantF) GetExpr(pos *int) Expr {
	if (*pos) == 0 {
		return c
	}
	(*pos)--
	return nil
}
func (c *ConstantF) SetExpr(pos *int, e Expr) (replace_me, replaced bool) {
	if (*pos) == 0 {
		return true, false
	}
	(*pos)--
	return false, false
}

func (s *System) GetExpr(pos *int) Expr {
	if (*pos) == 0 {
		return s
	}
	(*pos)--
	return nil
}
func (s *System) SetExpr(pos *int, e Expr) (replace_me, replaced bool) {
	if (*pos) == 0 {
		return true, false
	}
	(*pos)--
	return false, false
}

func (u *Neg) GetExpr(pos *int) Expr {
	if (*pos) == 0 {
		return u
	}
	(*pos)--
	return u.C.GetExpr(pos)
}
func (u *Neg) SetExpr(pos *int, e Expr) (replace_me, replaced bool) {
	if (*pos) == 0 {
		return true, false
	}
	(*pos)--
	if *pos == 0 {
		u.C = e
		return false, true
	}
	rme, repd := u.C.SetExpr(pos, e)
	if repd {
		return false, true
	}
	if rme {
		u.C = e
		return false, true
	}
	return false, repd
}

func (u *Abs) GetExpr(pos *int) Expr {
	if (*pos) == 0 {
		return u
	}
	(*pos)--
	return u.C.GetExpr(pos)
}
func (u *Abs) SetExpr(pos *int, e Expr) (replace_me, replaced bool) {
	if (*pos) == 0 {
		return true, false
	}
	(*pos)--
	if *pos == 0 {
		u.C = e
		return false, true
	}
	rme, repd := u.C.SetExpr(pos, e)
	if repd {
		return false, true
	}
	if rme {
		u.C = e
		return false, true
	}
	return false, repd
}

func (u *Sqrt) GetExpr(pos *int) Expr {
	if (*pos) == 0 {
		return u
	}
	(*pos)--
	return u.C.GetExpr(pos)
}
func (u *Sqrt) SetExpr(pos *int, e Expr) (replace_me, replaced bool) {
	if (*pos) == 0 {
		return true, false
	}
	(*pos)--
	if *pos == 0 {
		u.C = e
		return false, true
	}
	rme, repd := u.C.SetExpr(pos, e)
	if repd {
		return false, true
	}
	if rme {
		u.C = e
		return false, true
	}
	return false, repd
}

func (u *Sin) GetExpr(pos *int) Expr {
	if (*pos) == 0 {
		return u
	}
	(*pos)--
	return u.C.GetExpr(pos)
}
func (u *Sin) SetExpr(pos *int, e Expr) (replace_me, replaced bool) {
	if (*pos) == 0 {
		return true, false
	}
	(*pos)--
	if *pos == 0 {
		u.C = e
		return false, true
	}
	rme, repd := u.C.SetExpr(pos, e)
	if repd {
		return false, true
	}
	if rme {
		u.C = e
		return false, true
	}
	return false, repd
}

func (u *Cos) GetExpr(pos *int) Expr {
	if (*pos) == 0 {
		return u
	}
	(*pos)--
	return u.C.GetExpr(pos)
}
func (u *Cos) SetExpr(pos *int, e Expr) (replace_me, replaced bool) {
	if (*pos) == 0 {
		return true, false
	}
	(*pos)--
	if *pos == 0 {
		u.C = e
		return false, true
	}
	rme, repd := u.C.SetExpr(pos, e)
	if repd {
		return false, true
	}
	if rme {
		u.C = e
		return false, true
	}
	return false, repd
}

func (u *Tan) GetExpr(pos *int) Expr {
	if (*pos) == 0 {
		return u
	}
	(*pos)--
	return u.C.GetExpr(pos)
}
func (u *Tan) SetExpr(pos *int, e Expr) (replace_me, replaced bool) {
	if (*pos) == 0 {
		return true, false
	}
	(*pos)--
	if *pos == 0 {
		u.C = e
		return false, true
	}
	rme, repd := u.C.SetExpr(pos, e)
	if repd {
		return false, true
	}
	if rme {
		u.C = e
		return false, true
	}
	return false, repd
}

func (u *Exp) GetExpr(pos *int) Expr {
	if (*pos) == 0 {
		return u
	}
	(*pos)--
	return u.C.GetExpr(pos)
}

func (u *Log) GetExpr(pos *int) Expr {
	if (*pos) == 0 {
		return u
	}
	(*pos)--
	return u.C.GetExpr(pos)
}

func (u *PowI) GetExpr(pos *int) Expr {
	if (*pos) == 0 {
		return u
	}
	(*pos)--
	return u.Base.GetExpr(pos)
}

func (u *PowF) GetExpr(pos *int) Expr {
	if (*pos) == 0 {
		return u
	}
	(*pos)--
	return u.Base.GetExpr(pos)
}

func (n *PowE) GetExpr(pos *int) Expr {
	if (*pos) == 0 {
		return n
	}
	(*pos)--
	tmp := n.Base.GetExpr(pos)
	if tmp != nil {
		return tmp
	}
	return n.Power.GetExpr(pos)
}

func (n *Div) GetExpr(pos *int) Expr {
	if (*pos) == 0 {
		return n
	}
	(*pos)--
	tmp := n.Numer.GetExpr(pos)
	if tmp != nil {
		return tmp
	}
	return n.Denom.GetExpr(pos)
}

func (n *Add) GetExpr(pos *int) Expr {
	if (*pos) == 0 {
		return n
	}
	(*pos)--
	for _, C := range n.CS {
		if C == nil {
			continue
		}
		tmp := C.GetExpr(pos)
		if tmp != nil {
			return tmp
		}
		if *pos < 0 {
			return nil
		}
	}
	return nil
}

func (n *Mul) GetExpr(pos *int) Expr {
	if (*pos) == 0 {
		return n
	}
	(*pos)--
	for _, C := range n.CS {
		if C == nil {
			continue
		}
		tmp := C.GetExpr(pos)
		if tmp != nil {
			return tmp
		}
		if *pos < 0 {
			return nil
		}
	}
	return nil
}

func SwapExpr(orig, newt Expr, pos int) (ret Expr) {
	//   fmt.Printf( "SWAP orig  %v\n", orig )
	p := pos
	//   oldt := orig.GetExpr(&p)
	//   fmt.Printf( "SWAP (%d)\n%v\n%v\n", pos, oldt, newt )
	rme, _ := orig.SetExpr(&p, newt)
	if rme {
		ret = newt
	}

	//   fmt.Printf( "SWAP ret  %v\n", ret )
	return
}

func (u *Exp) SetExpr(pos *int, e Expr) (replace_me, replaced bool) {
	if (*pos) == 0 {
		return true, false
	}
	(*pos)--
	if *pos == 0 {
		u.C = e
		return false, true
	}
	rme, repd := u.C.SetExpr(pos, e)
	if repd {
		return false, true
	}
	if rme {
		u.C = e
		return false, true
	}
	return false, repd
}
func (u *Log) SetExpr(pos *int, e Expr) (replace_me, replaced bool) {
	if (*pos) == 0 {
		return true, false
	}
	(*pos)--
	if *pos == 0 {
		u.C = e
		return false, true
	}
	rme, repd := u.C.SetExpr(pos, e)
	if repd {
		return false, true
	}
	if rme {
		u.C = e
		return false, true
	}
	return false, repd
}
func (u *PowI) SetExpr(pos *int, e Expr) (replace_me, replaced bool) {
	if (*pos) == 0 {
		return true, false
	}
	(*pos)--
	if *pos == 0 {
		u.Base = e
		return false, true
	}
	rme, repd := u.Base.SetExpr(pos, e)
	if repd {
		return false, true
	}
	if rme {
		u.Base = e
		return false, true
	}
	return false, repd
}
func (u *PowF) SetExpr(pos *int, e Expr) (replace_me, replaced bool) {
	if (*pos) == 0 {
		return true, false
	}
	(*pos)--
	if *pos == 0 {
		u.Base = e
		return false, true
	}
	rme, repd := u.Base.SetExpr(pos, e)
	if repd {
		return false, true
	}
	if rme {
		u.Base = e
		return false, true
	}
	return false, repd
}

func (n *Add) SetExpr(pos *int, e Expr) (replace_me, replaced bool) {
	if (*pos) == 0 {
		return true, false
	}
	(*pos)--
	var rme, repd bool
	for i, C := range n.CS {
		if C == nil {
			continue
		}
		rme, repd = C.SetExpr(pos, e)
		if repd {
			return false, true
		}
		if rme {
			n.CS[i] = e
			return false, true
		}
	}
	if *pos == 0 {
		n.CS = append(n.CS, e)
		return false, true
	}
	return false, repd
}

func (n *Mul) SetExpr(pos *int, e Expr) (replace_me, replaced bool) {
	if (*pos) == 0 {
		return true, false
	}
	(*pos)--
	var rme, repd bool
	for i, C := range n.CS {
		if C == nil {
			continue
		}
		rme, repd = C.SetExpr(pos, e)
		if repd {
			return false, true
		}
		if rme {
			n.CS[i] = e
			return false, true
		}
	}
	if *pos == 0 {
		n.CS = append(n.CS, e)
		return false, true
	}
	return false, repd
}
func (n *Div) SetExpr(pos *int, e Expr) (replace_me, replaced bool) {
	if (*pos) == 0 {
		return true, false
	}
	(*pos)--
	if *pos == 0 {
		n.Numer = e
		return false, true
	}
	rme, repd := n.Numer.SetExpr(pos, e)
	if repd {
		return false, true
	}
	if rme {
		n.Numer = e
		return false, true
	}
	if *pos == 0 {
		n.Denom = e
		return false, true
	}
	rme, repd = n.Denom.SetExpr(pos, e)
	if repd {
		return false, true
	}
	if rme {
		n.Denom = e
		return false, true
	}
	return false, repd
}
func (n *PowE) SetExpr(pos *int, e Expr) (replace_me, replaced bool) {
	if (*pos) == 0 {
		return true, false
	}
	(*pos)--
	if *pos == 0 {
		n.Base = e
		return false, true
	}
	rme, repd := n.Base.SetExpr(pos, e)
	if repd {
		return false, true
	}
	if rme {
		n.Power = e
		return false, true
	}
	if *pos == 0 {
		n.Power = e
		return false, true
	}
	rme, repd = n.Power.SetExpr(pos, e)
	if repd {
		return false, true
	}
	if rme {
		n.Power = e
		return false, true
	}
	return false, repd
}
