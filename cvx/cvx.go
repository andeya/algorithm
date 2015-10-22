henrylee2cn/algorithm// Copyright (c) Harri Rautila, 2012

// This file is part of github.com/hrautila/cvx package. 
// It is free software, distributed under the terms of GNU Lesser General Public 
// License Version 3, or any later version. See the COPYING tile included in this archive.

package cvx

import (
    "github.com/henrylee2cn/algorithm/cvx/sets"
    "github.com/henrylee2cn/algorithm/matrix"
)

// kktFactor produces solver function
type kktFactor func(*sets.FloatMatrixSet, *matrix.FloatMatrix, *matrix.FloatMatrix) (KKTFunc, error)

// kktSolver creates problem spesific factor
type kktSolver func(*matrix.FloatMatrix, *sets.DimensionSet, *matrix.FloatMatrix, int) (kktFactor, error)

// Custom solver type for solving linear equations (`KKT systems')
//        
//            [ 0  A'  G'   ] [ ux ]   [ bx ]
//            [ A  0   0    ] [ uy ] = [ by ].
//            [ G  0  -W'*W ] [ uz ]   [ bz ]
//
// W is a scaling matrix, a block diagonal mapping 
//
// The call f = KKTConeSolver(W) or f = KKTCpSolver(W, x, z)should return a function f
// that solves the KKT system by f(x, y, z).  On entry, x, y, z contain the 
// righthand side bx, by, bz.  On exit, they contain the solution, 
// with uz scaled: the argument z contains W*uz.  In other words,
// on exit, x, y, z are the solution of
//
//            [ 0  A'  G'*W^{-1} ] [ ux ]   [ bx ]
//            [ A  0   0         ] [ uy ] = [ by ].
//            [ G  0  -W'        ] [ uz ]   [ bz ]
//

// KKTFunc solves KKT equations for matrix arguments
type KKTFunc func(x, y, z *matrix.FloatMatrix) error

// KKTFuncVar solves KKT equations for custom variable arguments
type KKTFuncVar func(x, y MatrixVariable, z *matrix.FloatMatrix) error

// KKTConeSolver produces solver function for cone problems with matrix variables
type KKTConeSolver func(W *sets.FloatMatrixSet) (KKTFunc, error)

// KKTConeSolver produces solver function for cone problems with custom variables
type KKTConeSolverVar func(W *sets.FloatMatrixSet) (KKTFuncVar, error)

// KKTCpSolver produces solver function for convex problems
type KKTCpSolver func(*sets.FloatMatrixSet, *matrix.FloatMatrix, *matrix.FloatMatrix) (KKTFunc, error)
type KKTCpSolverVar func(W *sets.FloatMatrixSet, x MatrixVariable, znl *matrix.FloatMatrix) (KKTFuncVar, error)

type solverMap map[string]kktSolver

var lpsolvers solverMap = solverMap{
    "ldl":   kktLdl,
    "ldl2":  kktLdl,
    "qr":    kktQr,
    "chol":  kktChol,
    "chol2": kktChol2}

var solvers solverMap = solverMap{
    "ldl":   kktLdl,
    "ldl2":  kktLdl,
    "chol":  kktChol,
    "chol2": kktChol2}

type StatusCode int

const (
    Optimal = StatusCode(1 + iota)
    PrimalInfeasible
    DualInfeasible
    Unknown
)

// If the exit status is Optimal, then the primal and dual
// infeasibilities are guaranteed to be less than 
// SolversOptions.FeasTol (default 1e-7).  The gap is less than
// SolversOptions.AbsTol (default 1e-7) or the relative gap is 
// less than SolversOptions.RelTol (defaults 1e-6).     
//
// Termination with status Unknown indicates that the algorithm 
// failed to find a solution that satisfies the specified tolerances.
// In some cases, the returned solution may be fairly accurate.  If
// the primal and dual infeasibilities, the gap, and the relative gap
// are small, then x, y, snl, sl, znl, zl are close to optimal.
//
type Solution struct {
    // Solution status
    Status StatusCode
    // Solution result set. 
    Result *sets.FloatMatrixSet
    // The primal objective c'*x
    PrimalObjective float64
    // The dual objective value
    DualObjective float64
    // Solution duality gap.
    Gap float64
    // Solution relative gap
    RelativeGap float64
    // Solution primal infeasibility
    PrimalInfeasibility float64
    // Solution dual infeasibility: the residual of the dual contraints
    DualInfeasibility float64
    // The smallest primal slack: min( min_k sl_k, sup{t | sl >= te}
    PrimalSlack float64
    // The smallest dual slack: min( min_k sl_k, sup{t | sl >= te}
    DualSlack          float64
    PrimalResidualCert float64
    DualResidualCert   float64
    // Number of iterations run
    Iterations int
}

// Solver options.
type SolverOptions struct {
    // Absolute tolerance
    AbsTol float64
    // Relative tolerance
    RelTol float64
    // Feasibility tolerance
    FeasTol float64
    // Maximum number of iterations
    MaxIter int
    // Show progress flag
    ShowProgress bool
    // Debug flag
    Debug bool
    // Refinement count
    Refinement int
    // KKT solver function name; "ldl", "ldl2", "qr", "chol", "chol2"
    // NOTE: currently all solvers mapped to "ldl".
    KKTSolverName string
}

const (
    MAXITERS = 100
    ABSTOL   = 1e-7
    RELTOL   = 1e-6
    FEASTOL  = 1e-7
)

// Local Variables:
// tab-width: 4
// End:
