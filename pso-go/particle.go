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
    "math/rand"
    "time"
)

// Range of values.
type Range struct {
	min []float64
	max []float64
}

// Create a new range.
func NewRange(min, max []float64) *Range {

	switch {
	case min == nil:
		panic("min cannot be nil.")
	case max == nil:
		panic("max cannot be nil.")
	case len(min) != len(max):
		panic("length of min and max have to be same.")
	}

	return &Range{min, max}
}

// Get either the vector is in this range or not.
func (r *Range) In(vector []float64) bool {

    switch {
    case vector == nil:
        panic("vector cannot be nil")
    }

	if len(vector) != len(r.min) {
		panic("length of values have to be same with minx and max.")
	}

	for i := range vector {
		if vector[i] < r.min[i] || vector[i] > r.max[i] {
			return false
		}
	}

	return true
}

func (r *Range) Min() []float64 {
    if r.min == nil {
        return nil
    }
    cpy := make([]float64, len(r.min))
    copy(cpy, r.min)
    return cpy
}

func (r *Range) Max() []float64 {
    if r.max == nil {
        return nil
    }
    cpy := make([]float64, len(r.max))
    copy(cpy, r.max)
    return cpy
}

// It is type of particle which find a result of optimization.
type Particle struct {
	// Current position
	position []float64
	// Current velocity
	velocity []float64
	// range of position
	valuesRange *Range
	// evaluated value by target function
	evalValue float64
	// local best of this particle
	best []float64
}

// Create a new particle.
func NewParticle(position, velocity []float64, valuesRange *Range) *Particle {

	switch {
	case position == nil:
		panic("position cannot be nil.")
	case velocity == nil:
		panic("velocity cannot be nil.")
    case len(position) != len(velocity):
        panic("length of position and velocity have to be same.")
	}

    cpyPos := make([]float64, len(position))
    copy(cpyPos, position)

    cpyVelocity := make([]float64, len(velocity))
    copy(cpyVelocity, velocity)


    best := make([]float64, len(position))
    copy(best, position)

	return &Particle{cpyPos, cpyVelocity, valuesRange, math.MaxFloat64, best}
}

// Get the position of particle on the solution space.
func (p *Particle) Position() []float64 {
	return p.position
}

// Get the position of particle.
func (p *Particle) Velocity() []float64 {
	return p.velocity
}

// Get the range of position.
func (p *Particle) Range() *Range {
	return p.valuesRange
}

// Get the evaluated value by target function.
func (p *Particle) EvalValue() float64 {
	return p.evalValue
}

// Get the local best of the particle.
func (p *Particle) Best() []float64 {
    if p.best == nil {
        return nil
    }

    cpy := make([]float64, len(p.best))
	return cpy
}

// Do a step of the particle.
func (p *Particle) Step(f TargetFunc, param *Param, globalBest []float64) {

    switch {
    case f == nil:
        panic("f cannot be nil.")
    case param == nil:
        panic("param cannot be nil.")
    case globalBest == nil:
        panic("globalBest cannot be nil.")
    case len(globalBest) != len(p.position):
        panic("length of particle position and globalBest have to be same")
    }

    oldPosition := make([]float64, len(p.position))
    c1 := param.C1()
    c2 := param.C2()
    w  := param.W()
    copy(oldPosition, p.position)

    // random generator
    var rnd = rand.New(rand.NewSource(time.Now().UnixNano()))

    for i := range p.position {
        // move
        p.position[i] += p.velocity[i]

        // Random value
        r1 := rnd.Float64()
        r2 := rnd.Float64()

        // update velocity
        p.velocity[i] = w[i] * p.velocity[i] + r1 * c1[i] * (p.best[i] - p.position[i]) * r2 * c2[i] * (globalBest[i] - p.position[i])
    }

    // Over the range?
    if !p.valuesRange.In(p.position) {
        copy(p.position, oldPosition)
    }

    // Update best
    p.evalValue = f(p.position)
    bestValue := f(p.best)
    if p.evalValue < bestValue {
        copy(p.best, p.position)
    }
}
