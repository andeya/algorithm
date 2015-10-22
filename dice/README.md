# Dice
## Simple dice rolling for Go.

`dice` offers three basic methods for rolling dice:
```go
dice.Roll()
```
,
```go
dice.RollP()
```    
, and
```go
dice.RollD6()
```
The first function receives a single `string` parameter which describes the dice roll.  For example, send `3d6` to roll three 6-sided dice, `1d20+1` to roll one 20-sided die and add 1 to the result, `1.5d6-1` to roll one and a half 6-sided rolls and subtract one from the result.  Half-dice are taken to be a single die of half the normal faces.  That is, `1.5d6-1` rolls one 6-sided die, one 3-sided die, and subtracts 1 from the sum of both dice.

The second function receives all the parameters of the dice roll separately (think of the `P` as "parameters", "pieces", or "parts").  The parameters are: `number` (the number of dice to roll), `faces` (how many sides, or faces, each die has), `adder` (any static quantity to add to the roll's total), and `half` (a boolean value indicating whether a half-die should be rolled in addition to the dice specified in the `number` parameter).  The equivalent to `dice.Roll("1.5d6-1")` is `dice.RollP(1,6,-1,true)`.

The third function receives a number and rolls that many 6-sided dice.  The last rolled die is considered a "wild" die.  A wild die coming up 6 causes an additional die to be rolled.  This additional die is also a wild die, which means that if it comes up 6, the process continues, until a non-6 die is rolled.

All three functions return a `DiceRoll` object:
```go
type DiceRoll struct {
    NumberOfDice int
    DieFaces     int
    Adder        int
    Half         bool
    Rolls        []int
    RawTotal     int
    Total        int
}
```    
`DiceRoll` objects have a method `Description()` which returns a description like the one used to roll with `dice.Roll()`.

Additionally, a `throwdice` runnable is included, which can take any number of dice roll descriptions as command-line arguments and prints the result for each. It will use the D6 system if passed the `-d6` flag.  It will also read from stdin if no command-line arguments are given, stopping on EOF or "exit" or "quit" (case insensitive).  Note that in D6 mode, descriptions are merely numbers, specifically the number of 6-sided dice to roll.

(Yes, it's all very simple and not extremely useful... it's mostly an exercise on Go on my part.)
