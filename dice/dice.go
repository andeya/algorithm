/*
Package dice provides simple dice-rolling methods allowing rolling of many varieties.

Three methods are exported: Roll(), RollP(), and RollD6().  Roll() receives a textual description
of the dice roll (e.g. "3d6+1"), while RollP() receives each parameter separately.
RollP() is faster, since it involves no parsing. RollD6() receives the number of 6-sided dice to
roll using the D6 system (see README for an explanation).
*/
package dice

import (
	"fmt"
	"regexp"
	"strconv"
)

/*
DiceRoll represents the result of a single dice roll.

    NumberOfDice       # how many dice are rolled
    DieFaces           # how many sides on each dice
    Adder              # static quantity to add to the total of the dice roll
    Half               # whether a half-die is rolled (in addition to
                             the number of dice specified in NumberOfDice)
    Rolls              # the individual dice rolls
    RawTotal           # the total of all the Rolls, the half-die (if any), and Adder
    Total              # The maximum of RawTotal and 1

Note that RawTotal can be negative (due to a negative Adder), but Total is never less than 1
*/
type DiceRoll struct {
	NumberOfDice int
	DieFaces     int
	Adder        int
	Half         bool
	Rolls        []int
	RawTotal     int
	Total        int
}

/*
Description returns a textual description of the dice roll.  It uses
the same format expected by Roll()
*/
func (dr *DiceRoll) Description() string {
	var adder, half string
	if dr.Adder < 0 {
		adder = strconv.Itoa(dr.Adder)
	} else if dr.Adder == 0 {
		adder = ""
	} else {
		adder = "+" + strconv.Itoa(dr.Adder)
	}
	if dr.Half {
		half = ".5"
	} else {
		half = ""
	}
	return fmt.Sprintf("%v%vd%v%v", dr.NumberOfDice, half, dr.DieFaces, adder)
}

func newDiceRollP(number, faces, adder int, half bool, rng intRng, repeatOnMaxLast bool) *DiceRoll {
	d := &DiceRoll{NumberOfDice: number, DieFaces: faces, Adder: adder, Half: half}
	d.Rolls = make([]int, number)
	for i := 0; i < number; i++ {
		d.Rolls[i] = rng.Intn(faces) + 1
		d.RawTotal += d.Rolls[i]
	}
	if repeatOnMaxLast {
		for d.Rolls[d.NumberOfDice-1] == faces {
			d.NumberOfDice = d.NumberOfDice + 1
			newDie := rng.Intn(faces) + 1
			d.Rolls = append(d.Rolls, newDie)
			d.RawTotal += newDie
		}
	}
	if half {
		halfDie := rng.Intn(faces/2) + 1
		d.Rolls = append(d.Rolls, halfDie)
		d.RawTotal += halfDie
	}
	d.RawTotal += d.Adder
	if d.RawTotal < 1 {
		d.Total = 1
	} else {
		d.Total = d.RawTotal
	}
	return d
}

/*
RollP() generates a new DiceRoll based on the specified parameters.
*/
func RollP(number, faces, adder int, half bool) *DiceRoll {
	return newDiceRollP(number, faces, adder, half, localRng, false)
}

/*
RollD6() generates a new DiceRoll based on the D6 system.

The number of dice specified is rolled, the last die being a "wild" die.  "Wild" dice
trigger an additional "wild" die to be rolled when coming up "6".
*/
func RollD6(number int) *DiceRoll {
	return newDiceRollP(number, 6, 0, false, localRng, true)
}

var diceExp = regexp.MustCompile(`^([1-9][0-9]*)?(\.5)?[dD]([1-9][0-9]*)([+-][1-9][0-9]*)?$`)

/*
Roll() generates a new dice roll based on parameters obtained from parsing the passed string.

The expected format is <NumDice>d<DieFaces>[Adder], where NumDice is the number of dice to roll,
DieFaces is how many sides per die, and Adder is a static quantity to add to the total.  Adder
must start with "+" or "-", and NumDice may end with ".5" to indicate a half die.
*/
func Roll(description string) (*DiceRoll, error) {
	return newDiceRoll(description, localRng, false)
}

func newDiceRoll(description string, rng intRng, repeatOnMaxLast bool) (*DiceRoll, error) {
	parts := diceExp.FindStringSubmatch(description)
	if parts == nil {
		return nil, fmt.Errorf("Bad description: %v", description)
	}

	numS, halfS, facesS, adderS := parts[1], parts[2], parts[3], parts[4]
	var number, faces, adder int
	var half bool
	var err error
	if numS == "" {
		number = 1
	} else if number, err = strconv.Atoi(numS); err != nil {
		return nil, fmt.Errorf("Bad number of dice: %v", numS)
	}
	half = halfS != ""
	if faces, err = strconv.Atoi(facesS); err != nil {
		return nil, fmt.Errorf("Bad die faces: %v", facesS)
	}
	if adderS == "" {
		adder = 0
	} else if adder, err = strconv.Atoi(adderS); err != nil {
		return nil, fmt.Errorf("Bad adder: %v", adderS)
	}
	return newDiceRollP(number, faces, adder, half, rng, repeatOnMaxLast), nil
}
