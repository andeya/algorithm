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
    "math"
)

// Target function
type TargetFunc func(vector []float64) float64

// Params of a solver.
type Param struct {
	w  []float64
	c1 []float64
	c2 []float64
}

// Create a new Param.
func NewParam(w, c1, c2 []float64) *Param {

	switch {
	case w == nil:
		panic("w cannot be nil.")
	case c1 == nil:
		panic("c1 cannot be nil.")
	case c2 == nil:
		panic("c2 cannot be nil.")
	case len(w) != len(c1) || len(w) != len(c2):
		panic("length of w and c1 and c2 have to be same.")
	}

	return &Param{w, c1, c2}
}

// Get w of param
func (p *Param) W() []float64 {
	return p.w
}

// Get c1 of param
func (p *Param) C1() []float64 {
	return p.c1
}

// Get c2 of param
func (p *Param) C2() []float64 {
	return p.c2
}

// Solver of PSO
type Solver struct {
	f      TargetFunc
	particles []*Particle
	param *Param
	best []float64
}

// Create a new solver.
func NewSolver(f TargetFunc, particles []*Particle, param *Param) *Solver {

	if particles == nil {
		panic("particles cannot be nil.")
	} else if len(particles) <= 0 {
		panic("The number of particles have to take more than 0.")
	}

    var bestValue float64
    var best []float64
    for _, p := range particles {
        if p == nil {
            continue
        }

        if best == nil {
            best = p.position
            bestValue = f(best)
        } else if f(p.position) < bestValue {
            copy(best, p.position)
            bestValue = f(best)
        }
    }

	return &Solver{f, particles, param, best}
}

// Get the target function.
func (s *Solver) TargetFunc() TargetFunc {
	return s.f
}

// Get the particles.
// It is a particles array which is copy of original one.
func (s *Solver) Particles() []*Particle {

	cpy := make([]*Particle, len(s.particles))
	copy(cpy, s.particles)

	return s.particles
}

// Get the parameter of the solver.
func (s *Solver) Param() *Param {
	return s.param
}

// Get best values.
// Before start steps it takes nil.
func (s *Solver) Best() []float64 {
    if s.best == nil {
        return nil
    }

    cpy := make([]float64, len(s.best))
    copy(cpy, s.best)
	return cpy
}

// Do a step of solver process.
func (s *Solver) Step() {

    bestValue := s.f(s.best)
    for _, p := range s.particles {
		if p == nil {
			continue
		}

        if s.best == nil {
            s.best = p.Best()
            bestValue = s.f(s.best)
        }

		p.Step(s.f, s.param, s.best)
        if p.EvalValue() < bestValue {
//            println("update!")
            copy(s.best, p.best)
            bestValue = p.EvalValue()
        }
	}

}

// Run solver process.
// If evaluated value is not change with error value or step count is over maxcount,
// the process will be stop.
func (s *Solver) Run(errorValue float64, maxcount uint) {
    pre := math.MaxFloat64
    for i := uint(0); i < maxcount; i++ {
        s.Step()
        v := s.f(s.best)
        if math.Abs(pre - v) <= errorValue {
            break
        }
        pre = v
    }
}
