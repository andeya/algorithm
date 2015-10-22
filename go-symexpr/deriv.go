package symexpr



func (l *Leaf) DerivVar( i int) Expr     { return NewConstantF(0.0) }
func (l *Leaf) DerivConst( i int) Expr  { return NewConstantF(0.0) }

func (u *Unary) DerivVar( i int ) Expr { return nil }
func (u *Unary) DerivConst( i int ) Expr { return nil }

func (n *N_ary) DerivVar( i int ) Expr { return nil }
func (n *N_ary) DerivConst( i int ) Expr { return nil }


func (v *Var) DerivVar( i int) Expr      {
  var f float64
  if v.P == i { f = 1.0 }
  return NewConstantF(f)
}


func (u *Neg) DerivVar( i int) Expr    { return NewNeg(u.C.DerivVar(i)) }
func (u *Abs) DerivVar( i int) Expr    { return NewAbs(u.C.DerivVar(i)) }
func (u *Sqrt) DerivVar( i int) Expr     { return NewPowF(u.C.Clone(),0.5).DerivVar(i) }
func (u *Sin) DerivVar( i int) Expr    {
  if u.C.HasVarI(i) {
    c := NewCos(u.C.Clone())
    g := u.C.DerivVar(i)
    m := NewMul()
    m.Insert(g)
    m.Insert(c)
    return m
  }
  return NewConstantF(0.0)
}
func (u *Cos) DerivVar( i int) Expr    {
  if u.C.HasVarI(i) {
    n := NewNeg(NewSin(u.C.Clone()))
    g := u.C.DerivVar(i)
    m := NewMul()
    m.Insert(g)
    m.Insert(n)
    return m
  }
  return NewConstantF(0.0)
}
func (u *Tan) DerivVar( i int) Expr    {
  if u.C.HasVarI(i) {
    n := u.C.DerivVar(i)
    d := NewPowI(NewCos(u.C.Clone()),2)
    return NewDiv(n,d)
  }
  return NewConstantF(0.0)
}
func (u *Exp) DerivVar( i int) Expr    {
  if u.C.HasVarI(i) {
    e := u.Clone()
    g := u.C.DerivVar(i)
    m := NewMul()
    m.Insert(g)
    m.Insert(e)
    return m
  }
  return NewConstantF(0.0)
}
func (u *Log) DerivVar( i int) Expr    {
  if u.C.HasVarI(i) {
    var d Div
    d.Numer = u.C.DerivVar(i)
    d.Denom = u.C.Clone()
    return &d
  }
  return NewConstantF(0.0)
}
func (u *PowI) DerivVar( i int) Expr     {
  if u.Base.HasVarI(i) {
    p := NewPowI(u.Base.Clone(),u.Power-1.0)
    c := &ConstantF{F: float64(u.Power)}
    g := u.Base.DerivVar(i)
    m := NewMul()
    m.Insert(c)
    m.Insert(g)
    m.Insert(p)
    return m
  }
  return NewConstantF(0.0)
}
func (u *PowF) DerivVar( i int) Expr     {
  if u.Base.HasVarI(i) {
    p := NewPowF(u.Base.Clone(),u.Power-1.0)
    c := &ConstantF{F: u.Power}
    g := u.Base.DerivVar(i)
    m := NewMul()
    m.Insert(c)
    m.Insert(g)
    m.Insert(p)
    return m
  }
  return NewConstantF(0.0)
}

func (n *Add) DerivVar( i int) Expr      {
  if n.HasVarI(i) {
    a := NewAdd()
    for _,C := range n.CS {
      if C == nil { continue }
      if C.HasVarI(i) {
        a.Insert( C.DerivVar(i) )
      }
    }
    if len(a.CS) > 0 {
      return a
    }
  }
  return NewConstantF(0.0)
}
func (n *Mul) DerivVar( i int) Expr      {
  if n.HasVarI(i) {
    a := NewAdd()
    for j,J := range n.CS {
      if J == nil { continue }
      if J.HasVarI(i) {
        m := NewMul()
        for I,C := range n.CS {
          if C == nil { continue }
//           fmt.Printf( "%d,%d  %v\n", j,I, C)
          if j==I {
            m.Insert( C.DerivVar(i) )
          } else {
            m.Insert( C.Clone() )
          }
        }
        a.Insert(m)
      }

    }
    if len(a.CS) > 0 {
      return a
    }
  }
  return NewConstantF(0.0)
}
func (n *Div) DerivVar( i int) Expr      {
  if n.HasVarI(i) {
    d := new(Div)

    a := NewAdd()
    m1 := NewMul()
    m1.Insert(n.Numer.DerivVar(i))
    m1.Insert(n.Denom.Clone())
    m2 := NewMul()
    m2.Insert(NewConstantF(-1.0))
    m2.Insert(n.Numer.Clone())
    m2.Insert(n.Denom.DerivVar(i))
    a.Insert(m1)
    a.Insert(m2)
    d.Numer = a

    p2 := NewPowI(n.Denom.Clone(),2)
    d.Denom = p2
    return d
  }
  return NewConstantF(0.0)
}

// TODO TODO TODO fix me and the next one
func (n *PowE) DerivVar( i int) Expr     {
  if n.HasVarI(i) {
    a := NewAdd()

    m1 := NewMul()
    m1.Insert(n.Base.DerivVar(i))
    m1.Insert(n.Power.Clone())

    // ???
    p1 := new(PowE)
    p1.Base = n.Base.Clone()
    // ???
    a1 := NewAdd()
    a1.Insert(n.Power.Clone())
    a1.Insert(NewConstantF(-1.0))

    m2 := NewMul()
    m2.Insert(n.Power.DerivVar(i))
    m2.Insert(n.Clone())
    m2.Insert(NewLog(n.Base.Clone()))

    a.Insert(m1)
    a.Insert(m2)


  }
  return NewConstantF(0.0)
}






func (c *Constant) DerivConst( i int) Expr {
  var f float64
  if c.P == i { f = 1.0 }
  return NewConstantF(f)
}


func (u *Neg) DerivConst( i int) Expr    { return NewNeg(u.C.DerivConst(i)) }
func (u *Abs) DerivConst( i int) Expr    { return NewAbs(u.C.DerivConst(i)) }
func (u *Sqrt) DerivConst( i int) Expr     { return (NewPowF(u.C.Clone(),0.5)).DerivConst(i) }
func (u *Sin) DerivConst( i int) Expr    {
  if u.C.HasConstI(i) {
    c := NewCos(u.C.Clone())
    g := u.C.DerivConst(i)
    m := NewMul()
    m.Insert(g)
    m.Insert(c)
    return m
  }
  return NewConstantF(0.0)
}
func (u *Cos) DerivConst( i int) Expr    {
  if u.C.HasConstI(i) {
    n := NewNeg(NewSin(u.C.Clone()))
    g := u.C.DerivConst(i)
    m := NewMul()
    m.Insert(g)
    m.Insert(n)
    return m
  }
  return NewConstantF(0.0)
}
func (u *Exp) DerivConst( i int) Expr    {
  if u.C.HasConstI(i) {
    e := u.Clone()
    g := u.C.DerivConst(i)
    m := NewMul()
    m.Insert(g)
    m.Insert(e)
    return m
  }
  return NewConstantF(0.0)
}
func (u *Log) DerivConst( i int) Expr    {
  if u.C.HasConstI(i) {
    return NewDiv(u.C.DerivConst(i),u.C.Clone())
  }
  return NewConstantF(0.0)
}
func (u *PowI) DerivConst( i int) Expr     {
  if u.Base.HasConstI(i) {
    p := NewPowI(u.Base.Clone(),u.Power-1)
    c := &ConstantF{F: float64(u.Power)}
    g := u.Base.DerivConst(i)
    m := NewMul()
    m.Insert(c)
    m.Insert(g)
    m.Insert(p)
    return m
  }
  return NewConstantF(0.0)
}
func (u *PowF) DerivConst( i int) Expr     {
  if u.Base.HasVarI(i) {
    p := NewPowF(u.Base.Clone(),u.Power-1.0)
    c := NewConstantF(u.Power)
    g := u.Base.DerivConst(i)
    m := NewMul()
    m.Insert(c)
    m.Insert(g)
    m.Insert(p)
    return m
  }
  return NewConstantF(0.0)
}

func (n *Add) DerivConst( i int) Expr      {
  if n.HasConstI(i) {
    a := NewAdd()
    for _,C := range n.CS {
      if C == nil { continue }
      if C.HasConstI(i) {
        a.Insert( C.DerivConst(i) )
      }
    }
    if len(a.CS) > 0 {
      return a
    }
  }
  return NewConstantF(0.0)
}

func (n *Mul) DerivConst( i int) Expr      {
  if n.HasConstI(i) {
    a := NewAdd()
    for j,J := range n.CS {
      if J == nil { continue }
      if J.HasConstI(i) {
        m := NewMul()
        for I,C := range n.CS {
          if C == nil { continue }
          if j==I {
            m.Insert( C.DerivConst(i) )
          } else {
            m.Insert( C.Clone() )
          }
        }
        a.Insert(m)
      }

    }
    if len(a.CS) > 0 {
      return a
    }
  }
  return NewConstantF(0.0)
}
func (n *Div) DerivConst( i int) Expr      {
  if n.HasConstI(i) {
    d := new(Div)

    a := NewAdd()
    m1 := NewMul()
    m1.Insert(n.Denom.Clone())
    m1.Insert(n.Numer.DerivConst(i))
    a.Insert(m1)

    m2 := NewMul()
    m2.Insert(n.Numer.Clone())
    m2.Insert(n.Denom.DerivConst(i))
    n2 := NewNeg(m2)
    a.Insert(n2)

    d.Numer = a

    p2 := NewPowI(n.Denom.Clone(),2)
    d.Denom = p2
    return d
  }
  return NewConstantF(0.0)
}
func (n *PowE) DerivConst( i int) Expr     {
  if n.HasConstI(i) {
    a := NewAdd()

    m1 := NewMul()
    m1.Insert(n.Base.DerivConst(i))
    m1.Insert(n.Power.Clone())
    p1 := new(PowE)
    p1.Base = n.Base.Clone()
    a1 := NewAdd()
    a1.Insert(n.Power.Clone())
    a1.Insert( NewConstantF(-1.0))

    m2 := NewMul()
    m2.Insert(n.Power.DerivConst(i))
    m2.Insert(n.Clone())
    m2.Insert( NewLog(n.Base.Clone()) )

    a.Insert(m1)
    a.Insert(m2)


  }
  return NewConstantF(0.0)
}
