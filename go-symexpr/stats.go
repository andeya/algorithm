package symexpr

type ExprParams struct {

	// bounds on tree
	MaxS, MaxD,
	MinS, MinD int

	// bounds on some operators
	NumDim, NumSys, NumCoeff int

	// usable terms at each location
	Roots, Nodes, Leafs, NonTrig []Expr
}

type ExprStats struct {
	depth   int // layers from top (root == 1)
	height  int // layers of subtree (leaf == 1)
	size    int
	pos     int
	numchld int
}

func (es *ExprStats) Depth() int       { return es.depth }
func (es *ExprStats) Height() int      { return es.height }
func (es *ExprStats) Size() int        { return es.size }
func (es *ExprStats) Pos() int         { return es.pos }
func (es *ExprStats) NumChildren() int { return es.numchld }

func (t *Leaf) CalcExprStats() {
	t.depth = 1
	t.height = 1
	t.size = 1
	t.pos = 0
	t.numchld = 0
}
func (t *Leaf) calcExprStatsR(depth int, pos *int) {
	t.depth = depth + 1
	t.height = 1
	t.size = 1
	t.pos = *pos
	(*pos)++
	t.numchld = 0
}

func (u *Unary) CalcExprStats() {
	pos := 1
	u.depth = 1
	u.C.calcExprStatsR(1, &pos)
	u.height = 1 + u.C.Height()
	u.size = 1 + u.C.Size()
	u.pos = 0
	u.numchld = 1
}
func (u *Unary) calcExprStatsR(depth int, pos *int) {
	u.pos = *pos
	(*pos)++
	u.depth = depth + 1
	u.C.calcExprStatsR(depth+1, pos)
	u.height = 1 + u.C.Height()
	u.size = 1 + u.C.Size()
	u.numchld = 1
}

func (n *N_ary) CalcExprStats() {
	c_high := 0
	pos := 1
	n.depth = 1
	n.size = 1
	n.pos = 0
	n.numchld = len(n.CS)
	for _, C := range n.CS {
		if C == nil {
			continue
		} else {
			n.numchld++
		}
		C.calcExprStatsR(1, &pos)
		n.size += C.Size()
		if c_high < C.Height() {
			c_high = C.Height()
		}
	}
	n.height = 1 + c_high
}

func (n *N_ary) calcExprStatsR(depth int, pos *int) {
	c_high := 0
	n.depth = depth + 1
	n.size = 1
	n.pos = *pos
	(*pos)++
	n.numchld = len(n.CS)
	for _, C := range n.CS {
		if C == nil {
			continue
		} else {
			n.numchld++
		}
		C.calcExprStatsR(depth+1, pos)
		n.size += C.Size()
		if c_high < C.Height() {
			c_high = C.Height()
		}
	}
	n.height = 1 + c_high
}

func (u *PowI) CalcExprStats() {
	pos := 1
	u.Base.calcExprStatsR(1, &pos)
	u.depth = 1
	u.height = 1 + u.Base.Height()
	u.size = 1 + u.Base.Size()
	u.pos = 0
	u.numchld = 1
}

func (u *PowI) calcExprStatsR(depth int, pos *int) {
	u.pos = *pos
	(*pos)++
	u.Base.calcExprStatsR(depth+1, pos)
	u.depth = depth + 1
	u.height = 1 + u.Base.Height()
	u.size = 1 + u.Base.Size()
	u.numchld = 1
}

func (u *PowF) CalcExprStats() {
	pos := 1
	u.Base.calcExprStatsR(1, &pos)
	u.depth = 1
	u.height = 1 + u.Base.Height()
	u.size = 1 + u.Base.Size()
	u.pos = 0
	u.numchld = 1
}

func (u *PowF) calcExprStatsR(depth int, pos *int) {
	u.pos = *pos
	(*pos)++
	u.Base.calcExprStatsR(depth+1, pos)
	u.depth = depth + 1
	u.height = 1 + u.Base.Height()
	u.size = 1 + u.Base.Size()
	u.numchld = 1
}

func max(l, r int) int {
	if l > r {
		return l
	}
	return r
}

func (n *PowE) CalcExprStats() {
	pos := 1
	n.Base.calcExprStatsR(1, &pos)
	n.Power.calcExprStatsR(1, &pos)
	n.depth = 1
	n.height = 1 + max(n.Base.Height(), n.Power.Height())
	n.size = 1 + n.Base.Size() + n.Power.Size()
	n.numchld = 2
}
func (n *PowE) calcExprStatsR(depth int, pos *int) {
	n.pos = *pos
	(*pos)++
	n.Base.calcExprStatsR(depth+1, pos)
	n.Power.calcExprStatsR(depth+1, pos)
	n.depth = depth + 1
	n.height = 1 + max(n.Base.Height(), n.Power.Height())
	n.size = 1 + n.Base.Size() + n.Power.Size()
	n.numchld = 2
}

func (n *Div) CalcExprStats() {
	pos := 1
	n.Numer.calcExprStatsR(1, &pos)
	n.Denom.calcExprStatsR(1, &pos)
	n.depth = 1
	n.height = 1 + max(n.Numer.Height(), n.Denom.Height())
	n.size = 1 + n.Numer.Size() + n.Denom.Size()
	n.numchld = 2
}
func (n *Div) calcExprStatsR(depth int, pos *int) {
	n.pos = *pos
	(*pos)++
	n.Numer.calcExprStatsR(depth+1, pos)
	n.Denom.calcExprStatsR(depth+1, pos)
	n.depth = depth + 1
	n.height = 1 + max(n.Numer.Height(), n.Denom.Height())
	n.size = 1 + n.Numer.Size() + n.Denom.Size()
	n.numchld = 2
}
