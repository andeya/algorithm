// Package evaler implements a simple fp arithmetic expression evaluator.
//
// See README.md for documentation.

package evaler

import (
	"fmt"
	"math"
	"math/big"
	"regexp"
	"strconv"
	"strings"

	"github.com/henrylee2cn/algorithm/evaler/stack"
)

var whitespace_rx = regexp.MustCompile(`\s+`)

// Unary minus is appeared at the following positions.
//     * the beginning of an expression
//     * after an operator or '('
var unary_minus_rx = regexp.MustCompile(`((?:^|[-+*/<>(])\s*)-`)
var fp_rx = regexp.MustCompile(`(\d+(?:\.\d+)?)`) // simple fp number

// Operator '@' means unary minus
var operators = "-+**/<>@"

// prec returns the operator's precedence
func prec(op string) (result int) {
	if op == "-" || op == "+" {
		result = 1
	} else if op == "*" || op == "/" {
		result = 2
	} else if op == "**" {
		result = 3
	} else if op == "@" {
		result = 4
	}
	return
}

// opGTE returns true if op1's precedence is >= op2
func opGTE(op1, op2 string) bool {
	return prec(op1) >= prec(op2)
}

// isOperator returns true if token is an operator
func isOperator(token string) bool {
	return strings.Contains(operators, token)
}

// isOperand returns true if token is an operand
func isOperand(token string) bool {
	return fp_rx.MatchString(token)
}

// convert2postfix converts an infix expression to postfix
func convert2postfix(tokens []string) []string {
	var stack stack.Stack
	var result []string
	for _, token := range tokens {

		if isOperator(token) {

		OPERATOR:
			for {
				top, err := stack.Top()
				if err == nil && top != "(" {
					if opGTE(top.(string), token) {
						pop, _ := stack.Pop()
						result = append(result, pop.(string))
					} else {
						break OPERATOR
					}
				}
				break OPERATOR
			}
			stack.Push(token)

		} else if token == "(" {
			stack.Push(token)

		} else if token == ")" {
		PAREN:
			for {
				top, err := stack.Top()
				if err == nil && top != "(" {
					pop, _ := stack.Pop()
					result = append(result, pop.(string))
				} else {
					stack.Pop() // pop off "("
					break PAREN
				}
			}

		} else if isOperand(token) {
			result = append(result, token)
		}

	}

	for !stack.IsEmpty() {
		pop, _ := stack.Pop()
		result = append(result, pop.(string))
	}

	return result
}

// evaluatePostfix takes a postfix expression and evaluates it
func evaluatePostfix(postfix []string) (*big.Rat, error) {
	var stack stack.Stack
	result := new(big.Rat) // note: a new(big.Rat) has value "0/1" ie zero
	for _, token := range postfix {
		if isOperand(token) {
			bigrat := new(big.Rat)
			if _, err := fmt.Sscan(token, bigrat); err != nil {
				return nil, fmt.Errorf("unable to scan %s", token)
			}
			stack.Push(bigrat)
		} else if isOperator(token) {

			op2, err2 := stack.Pop()
			if err2 != nil {
				return nil, err2
			}

			var op1 interface{}
			if token != "@" {
				var err1 error
				if op1, err1 = stack.Pop(); err1 != nil {
					return nil, err1
				}
			}

			dummy := new(big.Rat)
			switch token {
			case "**":
				float1 := BigratToFloat(op1.(*big.Rat))
				float2 := BigratToFloat(op2.(*big.Rat))
				float_result := math.Pow(float1, float2)
				stack.Push(FloatToBigrat(float_result))
			case "*":
				result := dummy.Mul(op1.(*big.Rat), op2.(*big.Rat))
				stack.Push(result)
			case "/":
				result := dummy.Quo(op1.(*big.Rat), op2.(*big.Rat))
				stack.Push(result)
			case "+":
				result = dummy.Add(op1.(*big.Rat), op2.(*big.Rat))
				stack.Push(result)
			case "-":
				result = dummy.Sub(op1.(*big.Rat), op2.(*big.Rat))
				stack.Push(result)
			case "<":
				if op1.(*big.Rat).Cmp(op2.(*big.Rat)) <= -1 {
					stack.Push(big.NewRat(1, 1))
				} else {
					stack.Push(new(big.Rat))
				}
			case ">":
				if op1.(*big.Rat).Cmp(op2.(*big.Rat)) >= 1 {
					stack.Push(big.NewRat(1, 1))
				} else {
					stack.Push(new(big.Rat))
				}
			case "@":
				result := dummy.Mul(big.NewRat(-1, 1), op2.(*big.Rat))
				stack.Push(result)
			}
		} else {
			return nil, fmt.Errorf("unknown token %v", token)
		}
	}

	retval, err := stack.Pop()
	if err != nil {
		return nil, err
	}
	return retval.(*big.Rat), nil
}

// tokenise takes an expr string and converts it to a slice of tokens
//
// tokenise puts spaces around all non-numbers, removes leading and
// trailing spaces, then splits on spaces
//
func tokenise(expr string) []string {
	spaced := unary_minus_rx.ReplaceAllString(expr, "$1 @")
	spaced = fp_rx.ReplaceAllString(spaced, " ${1} ")
	symbols := []string{"(", ")"}
	for _, symbol := range symbols {
		spaced = strings.Replace(spaced, symbol, fmt.Sprintf(" %s ", symbol), -1)
	}
	stripped := whitespace_rx.ReplaceAllString(strings.TrimSpace(spaced), "|")
	return strings.Split(stripped, "|")
}

// Eval takes an infix string arithmetic expression, and evaluates it
//
// Usage:
//   result, err := evaler.Eval("1+2")
// Returns: the result of the evaluation, and any errors
//
func Eval(expr string) (result *big.Rat, err error) {
	defer func() {
		if e := recover(); e != nil {
			result = nil
			err = fmt.Errorf("Invalid Expression: %s", expr)
		}
	}()

	tokens := tokenise(expr)
	postfix := convert2postfix(tokens)
	return evaluatePostfix(postfix)
}

// BigratToInt converts a *big.Rat to an int64 (with truncation); it
// returns an error for integer overflows.
func BigratToInt(bigrat *big.Rat) (int64, error) {
	float_string := bigrat.FloatString(0)
	return strconv.ParseInt(float_string, 10, 64)
}

// BigratToInt converts a *big.Rat to a *big.Int (with truncation)
func BigratToBigint(bigrat *big.Rat) *big.Int {
	int_string := bigrat.FloatString(0)
	bigint := new(big.Int)
	// no error scenario could be imagined in testing, so discard err
	fmt.Sscan(int_string, bigint)
	return bigint
}

// BigratToFloat converts a *big.Rat to a float64 (with loss of
// precision).
func BigratToFloat(bigrat *big.Rat) float64 {
	float_string := bigrat.FloatString(10) // arbitrary largish precision
	// no error scenario could be imagined in testing, so discard err
	float, _ := strconv.ParseFloat(float_string, 64)
	return float
}

// FloatToBigrat converts a float64 to a *big.Rat.
func FloatToBigrat(float float64) *big.Rat {
	float_string := fmt.Sprintf("%g", float)
	bigrat := new(big.Rat)
	// no error scenario could be imagined in testing, so discard err
	fmt.Sscan(float_string, bigrat)
	return bigrat
}
