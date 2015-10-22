
Symbolic Expressions (just math for now)
========================================
 - represetation: as an AST tree
 - manipulation: add,rm, simplify, derivative
 - evaluation: explicit, RK4 integration, non-linear regression via levmar->MINPACK
 - printing: String() print, prettyprint, serialization

``` go
type exprType int

const (
	NULL      exprType = iota
	CONSTANT           // indexed constant, useful for regression tasks
	CONSTANTF          // floating point constant
	TIME               // useful when looking at time series and RK4 integration
	SYSTEM             // i use this like a variable that changes between experiments, but not with time (mass,size,etc.)
	VAR                // a canonical variable

	NEG
	ABS
	SQRT
	SIN
	COS
	TAN
	EXP
	LOG
	POWI // Expr^Integer
	POWF // Expr^Float

	POWE // Expr^Expr
	DIV

	ADD // these can have more than two child nodes
	MUL // this eases simplification

	EXPR_MAX
	STARTVAR // for serialization reduction of variables
)

// Expr is the interface to all node types for the AST of mathematical expression
// 

type Expr interface {

	// types.go (this file)
	ExprType() exprType
	Clone() Expr

	// stats.go
	Size() int
	Depth() int
	Height() int
	NumChildren() int
	CalcExprStats(currDepth int) (mySize int)

	// compare.go
	AmILess(rhs Expr) bool
	AmIEqual(rhs Expr) bool
	AmISame(rhs Expr) bool       // equality without coefficient values/index
	AmIAlmostSame(rhs Expr) bool // adds flexibility to mul comparison to AmISame
	Sort()

	// has.go
	HasVar() bool
	HasVarI(i int) bool
	NumVar() int

	// DFS for a floating point valued ConstantF
	HasConst() bool
	// DFS for a indexed valued Constant
	HasConstI(i int) bool
	// Counts the number of indexed Constant nodes
	NumConstants() int

	// convert.go

	// Converts indexed Constant nodes to ConstantF nodes
	// using the input slice as the values for replacement 
	ConvertToConstantFs(cs []float64) Expr
	// DFS converting float valued constants to indexed constants
	// the input should be an empty slice
	// the output is an appended slice the size = |ConstantF| in the tree
	ConvertToConstants(cs []float64) ([]float64, Expr)
	//   IndexConstants( ci int ) int

	// getset.go   
	// DFS retrieval of a node by index
	GetExpr(pos *int) Expr
	// DFS replacement of a node and it's subtree
	// replaced is used to discontinue the DFS after replacement
	// replace_me gets triggered when pos == 0 and informs the parent node to replace the respective child node
	SetExpr(pos *int, e Expr) (replace_me, replaced bool)

	// print.go

	// prints the AST 
	String() string

	// creates an integer representation of the AST in ~prefix notation
	// The input is an empty slice, output is the representation.
	// The output is generally the ExprType integer value
	// Associative operators (+ & *) also include the number of children.
	// The terminal nodes include the index when appropriate.
	Serial([]int) []int

	// Pretty print acts like String, but replaces the internal indexed
	// formatting with user specified strings and values
	PrettyPrint(dnames, snames []string, cvals []float64) string
	// 	WriteString( buf *bytes.Buffer )

	// eval.go
	// Evaluates an expression at one point
	// t is a time value
	// x are the input Var values
	// c are the indexed Constant values
	// s are the indexed System values
	// the output is the result of DFS evaluation
	Eval(t float64, x, c, s []float64) float64

	// simp.go
	Simplify(rules SimpRules) Expr

	// deriv.go
	// Calculate the derivative w.r.t. Var_i
	DerivVar(i int) Expr
	// Calculate the derivative w.r.t. Constant_i
	DerivConst(i int) Expr
}
```

Tony Worm  Sept, 2012
