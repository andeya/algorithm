package dice

import (
	"fmt"
	"sync"
	"testing"
)

type testRng struct {
	Max     int
	current int
	lk      sync.Mutex
}

func (r *testRng) Intn(n int) (v int) {
	r.lk.Lock()
	v = r.current
	r.current = (r.current + 1) % r.Max
	r.lk.Unlock()
	if v > n-1 {
		v = n - 1
	}
	return
}

func TestRoll(t *testing.T) {
	rng := &testRng{Max: 6}
	dice := newDiceRollP(6, 6, -1, false, rng, false)
	if dice.Total != 20 {
		t.Errorf("Expected 20, got %v.", dice.Total)
	}
	rng = &testRng{Max: 10}
	dice = newDiceRollP(20, 10, 2, false, rng, false)
	if dice.Total != 112 {
		t.Errorf("Expected 112, got %v.", dice.Total)
	}
	rng = &testRng{Max: 6}
	dice = newDiceRollP(1, 6, 0, true, rng, false)
	if dice.Total != 3 {
		t.Errorf("Expected 3, got %v.", dice.Total)
	}
}

func TestRollD6(t *testing.T) {
	rng := &testRng{Max: 6}
	dice := newDiceRollP(6, 6, 0, false, rng, true)
	if dice.NumberOfDice != 7 {
		t.Errorf("Expected 7 dice, got %v.", dice.NumberOfDice)
	}
	if dice.Total != 22 {
		t.Errorf("Expected 22, got %v.", dice.Total)
	}
	rng = &testRng{Max: 6}
	dice = newDiceRollP(1, 6, 0, false, rng, true)
	if dice.NumberOfDice != 1 {
		t.Errorf("Expected 1 die, got %v.", dice.NumberOfDice)
	}
	if dice.Total != 1 {
		t.Errorf("Expected 1, got %v.", dice.Total)
	}
}

func ExampleDiceRoll_Description() {
	dice := RollP(3, 6, -1, false)
	fmt.Println(dice.Description())
	dice = RollP(2, 6, 0, true)
	fmt.Println(dice.Description())
	// Output: 3d6-1
	// 2.5d6
}

func ExampleDiceRoll_Description_text() {
	dice, _ := Roll("3d6-2")
	fmt.Println(dice.Description())
	fmt.Println(dice.NumberOfDice)
	fmt.Println(dice.DieFaces)
	fmt.Println(dice.Adder)
	// Output: 3d6-2
	// 3
	// 6
	// -2
}

func BenchmarkRoll(b *testing.B) {
	for i := 0; i < b.N; i++ {
		_, _ = Roll("3d6+1")
	}
}

func BenchmarkRollP(b *testing.B) {
	for i := 0; i < b.N; i++ {
		_ = RollP(3, 6, 1, false)
	}
}
