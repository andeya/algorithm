/*
pso-go - PSO (Particle Swarm Optimization) library for Go.
https://github.com/tenntenn/pso-go

Copyright (c) 2012, Takuya Ueda.
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice,
  this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.
* Neither the name of the author nor the names of its contributors may be used
  to endorse or promote products derived from this software
  without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

package pso

import (
	"testing"
    "math"
    "math/rand"
    "time"
)

type Pos []float64

func NewPos(x, y float64) Pos {
    return Pos([]float64{x, y})
}

func (p Pos) X() float64 {
    return p[0]
}

func (p Pos) Y() float64 {
    return p[1]
}

func TestSolvingSimultaneousEquation(t *testing.T) {

    t.Log("Test solving simultaneous equation as followings.")
    t.Log("x + y - 3 = 0")
    t.Log("2x + 5y - 9 = 0")

	// x + y - 3 = 0
	// 2x + 5y - 9 = 0
	f := func(vector []float64) float64 {
        p := Pos(vector)

		return math.Pow(p.X() + p.Y() - 3, 2) + math.Pow(2 * p.X() + 5 * p.Y() - 9, 2)
	}

    // random generator
    var rnd = rand.New(rand.NewSource(time.Now().UnixNano()))

	// Create particles
	const count = 1000
	particles := make([]*Particle, count)
	for i := range particles {
		position := NewPos(float64(rnd.Intn(20) - 10), float64(rnd.Intn(20) - 10))
		velocity := NewPos(float64(rnd.Intn(20) - 10), float64(rnd.Intn(20) - 10))
        min := NewPos(-10, -10)
        max := NewPos(10, 10)
		valuesRange := NewRange(min, max)
		particles[i] = NewParticle(position, velocity, valuesRange)
	}

	// Create a solver
	w := NewPos(0.9, 0.9)
	c1 := NewPos(0.9, 0.9)
	c2 := NewPos(0.9, 0.9)
	param := NewParam(w, c1, c2)
	solver := NewSolver(TargetFunc(f), particles, param)

    const STEP = 1000
    solver.Run(0.00001, STEP)

    best := Pos(solver.Best())

    // check X
    if math.Abs(best.X() - 2.0) >= 0.001 {
        t.Errorf("Expect x = 2.0 but accutual %f", best.X())
    }

    // check Y
    if math.Abs(best.Y() - 1.0) >= 0.001 {
        t.Errorf("Expect y = 1.0 but accutual %f", best.X())
    }
}
