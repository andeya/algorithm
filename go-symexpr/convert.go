package symexpr



func (n *Time) ConvertToConstantFs( cs []float64 ) Expr     { return n }
func (v *Var) ConvertToConstantFs( cs []float64 ) Expr      { return v }
func (c *Constant) ConvertToConstantFs( cs []float64 ) Expr { return &ConstantF{F: cs[c.P]} }
func (c *ConstantF) ConvertToConstantFs( cs []float64 ) Expr { return c }
func (s *System) ConvertToConstantFs( cs []float64 ) Expr   { return s }

func (u *Neg) ConvertToConstantFs( cs []float64 ) Expr     {
  e := u.C.ConvertToConstantFs(cs)
  if u.C != e {
    u.C = e
  }
  return u
}
func (u *Abs) ConvertToConstantFs( cs []float64 ) Expr     {
  e := u.C.ConvertToConstantFs(cs)
  if u.C != e {
    u.C = e
  }
  return u
}
func (u *Sqrt) ConvertToConstantFs( cs []float64 ) Expr    {
  e := u.C.ConvertToConstantFs(cs)
  if u.C != e {
    u.C = e
  }
  return u
}
func (u *Sin) ConvertToConstantFs( cs []float64 ) Expr     {
  e := u.C.ConvertToConstantFs(cs)
  if u.C != e {
    u.C = e
  }
  return u
}
func (u *Cos) ConvertToConstantFs( cs []float64 ) Expr     {
  e := u.C.ConvertToConstantFs(cs)
  if u.C != e {
    u.C = e
  }
  return u
}
func (u *Tan) ConvertToConstantFs( cs []float64 ) Expr     {
  e := u.C.ConvertToConstantFs(cs)
  if u.C != e {
    u.C = e
  }
  return u
}
func (u *Exp) ConvertToConstantFs( cs []float64 ) Expr     {
  e := u.C.ConvertToConstantFs(cs)
  if u.C != e {
    u.C = e
  }
  return u
}
func (u *Log) ConvertToConstantFs( cs []float64 ) Expr     {
  e := u.C.ConvertToConstantFs(cs)
  if u.C != e {
    u.C = e
  }
  return u
}
func (u *PowI) ConvertToConstantFs( cs []float64 ) Expr    {
  e := u.Base.ConvertToConstantFs(cs)
  if u.Base != e {
    u.Base = e
  }
  return u
}
func (u *PowF) ConvertToConstantFs( cs []float64 ) Expr    {
  e := u.Base.ConvertToConstantFs(cs)
  if u.Base != e {
    u.Base = e
  }
  return u
}

func (n *Add) ConvertToConstantFs( cs []float64 ) Expr      {
  for i,_ := range n.CS {
    if n.CS[i] != nil {
      e := n.CS[i].ConvertToConstantFs(cs)
      if n.CS[i] != e {
        n.CS[i] = e
      }
    }
  }
  return n
}

func (n *Mul) ConvertToConstantFs( cs []float64 ) Expr      {
  for i,_ := range n.CS {
    if n.CS[i] != nil {
      e := n.CS[i].ConvertToConstantFs(cs)
      if n.CS[i] != e {
        n.CS[i] = e
      }
    }
  }
  return n
}

func (n *Div) ConvertToConstantFs( cs []float64 ) Expr      {
  e1,e2 := n.Numer.ConvertToConstantFs(cs), n.Denom.ConvertToConstantFs(cs)
  if n.Numer != e1 { n.Numer = e1 }
  if n.Denom != e2 { n.Denom = e2 }
  return n

}
func (n *PowE) ConvertToConstantFs( cs []float64 ) Expr     {
  e1,e2 := n.Base.ConvertToConstantFs(cs), n.Power.ConvertToConstantFs(cs)
  if n.Base != e1 { n.Base = e1 }
  if n.Power != e2 { n.Power = e2 }
  return n
}






func (n *Time) ConvertToConstants( cs []float64 ) ( []float64, Expr )     { return cs,n }
func (v *Var) ConvertToConstants( cs []float64 ) ( []float64, Expr )      { return cs,v }
func (c *Constant) ConvertToConstants( cs []float64 ) ( []float64, Expr ) {
  c.P = len(cs)
  return append(cs,float64(c.P)),c
}
func (c *ConstantF) ConvertToConstants( cs []float64 ) ( []float64, Expr ) {
  C := &Constant{P:len(cs)}
  return append(cs,c.F),C
}
func (s *System) ConvertToConstants( cs []float64 ) ( []float64, Expr )   { return cs,s }

func (u *Neg) ConvertToConstants( cs []float64 ) ( []float64, Expr )     {
  css,e := u.C.ConvertToConstants(cs)
  if u.C != e { u.C = e }
  return css,u
}
func (u *Abs) ConvertToConstants( cs []float64 ) ( []float64, Expr )     {
  css,e := u.C.ConvertToConstants(cs)
  if u.C != e { u.C = e }
  return css,u
}
func (u *Sqrt) ConvertToConstants( cs []float64 ) ( []float64, Expr )    {
  css,e := u.C.ConvertToConstants(cs)
  if u.C != e { u.C = e }
  return css,u
}
func (u *Sin) ConvertToConstants( cs []float64 ) ( []float64, Expr )     {
  css,e := u.C.ConvertToConstants(cs)
  if u.C != e { u.C = e }
  return css,u
}
func (u *Cos) ConvertToConstants( cs []float64 ) ( []float64, Expr )     {
  css,e := u.C.ConvertToConstants(cs)
  if u.C != e { u.C = e }
  return css,u
}
func (u *Tan) ConvertToConstants( cs []float64 ) ( []float64, Expr )     {
  css,e := u.C.ConvertToConstants(cs)
  if u.C != e { u.C = e }
  return css,u
}
func (u *Exp) ConvertToConstants( cs []float64 ) ( []float64, Expr )     {
  css,e := u.C.ConvertToConstants(cs)
  if u.C != e { u.C = e }
  return css,u
}
func (u *Log) ConvertToConstants( cs []float64 ) ( []float64, Expr )     {
  css,e := u.C.ConvertToConstants(cs)
  if u.C != e { u.C = e }
  return css,u
}
func (u *PowI) ConvertToConstants( cs []float64 ) ( []float64, Expr )    {
  css,e := u.Base.ConvertToConstants(cs)
  if u.Base != e { u.Base = e }
  return css,u
}
func (u *PowF) ConvertToConstants( cs []float64 ) ( []float64, Expr )    {
  css,e := u.Base.ConvertToConstants(cs)
  if u.Base != e { u.Base = e }
  return css,u
}

func (n *Add) ConvertToConstants( cs []float64 ) ( []float64, Expr )      {
  for i,_ := range n.CS {
    if n.CS[i] != nil {
      var e Expr
      cs,e = n.CS[i].ConvertToConstants(cs)
      if n.CS[i] != e { n.CS[i] = e }
    }
  }
  return cs,n
}


func (n *Mul) ConvertToConstants( cs []float64 ) ( []float64, Expr )      {
  for i,_ := range n.CS {
    if n.CS[i] != nil {
      var e Expr
      cs,e = n.CS[i].ConvertToConstants(cs)
      if n.CS[i] != e { n.CS[i] = e }
    }
  }
  return cs,n
}

func (n *Div) ConvertToConstants( cs []float64 ) ( []float64, Expr )      {
  var e Expr
  cs,e = n.Numer.ConvertToConstants(cs)
  if n.Numer != e { n.Numer = e }
  cs,e = n.Denom.ConvertToConstants(cs)
  if n.Denom != e { n.Denom = e }
  return cs,n
}
func (n *PowE) ConvertToConstants( cs []float64 ) ( []float64, Expr )     {
  var e Expr
  cs,e = n.Base.ConvertToConstants(cs)
  if n.Base != e { n.Base = e }
  cs,e = n.Power.ConvertToConstants(cs)
  if n.Power != e { n.Power = e }
  return cs,n
}
