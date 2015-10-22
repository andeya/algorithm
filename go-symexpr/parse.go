package symexpr

import (
	"fmt"
	"math"
	"strconv"
	"strings"
)

func ParseFunc(text string, varNames []string) Expr {
	expr := parse(text, varNames)
	return expr
}

var itemPrec = map[itemType]int{
	itemAdd:     3,
	itemNeg:     4,
	itemMul:     5,
	itemDiv:     5,
	itemDiv2:    5,
	itemKeyword: 6,
	itemCarrot:  7,
}

func parse(input string, variables []string) Expr {
	L := lex("test", input, variables)
	go L.run()

	// var rules = DefaultRules()

	// fmt.Printf( "Input:     %s\n",input)
	expr := parseExpr("", L, 0)
	// fmt.Printf( "Result:    %v\n\n",expr)
	// simp := expr.Simplify(rules)
	// fmt.Printf( "Output:    %v\n\n",simp)

	return expr
}

func parseExpr(prefix string, L *lexer, p int) Expr {
	// fmt.Printf("%sin-E(%d)\n", prefix, p)
	e := parsePiece(prefix+"  ", L)
	// fmt.Printf("%se(%d): %v\n", prefix, p, e)
	for i := 0; ; i++ {
		var next item
		// fmt.Printf("%sLEX: %v\n", prefix, L.nextToken)

		if L.nextToken.typ == itemNil {
			next = <-L.items
		} else {
			next, L.nextToken.typ = L.nextToken, itemNil
		}
		// fmt.Printf("%snext(%d-%d): %v\n", prefix, p, i, next)
		typ := next.typ
		// fmt.Printf("%sn(%d): %v\n", prefix, p, next)
		if isBinary(typ) && itemPrec[typ] >= p {
			// fmt.Printf("%sbin: %v\n", prefix, next)
			q := itemPrec[typ]
			e2 := parseExpr(prefix+"  ", L, q)
			// fmt.Printf("%se2: %v\n", prefix, e2)
			// fmt.Printf("%snext2(%d-%d): %v\n", prefix, p, i, next2)
			switch typ {
			case itemAdd:
				add := NewAdd()
				add.Insert(e)
				add.Insert(e2)
				e = add
			case itemMul:
				mul := NewMul()
				mul.Insert(e)
				mul.Insert(e2)
				e = mul
			case itemNeg:
				if e == nil {
					neg := NewNeg(e2)
					e = neg
				} else {
					add := NewAdd()
					add.Insert(e)
					add.Insert(NewNeg(e2))
					e = add
				}
			case itemDiv, itemDiv2:
				div := NewDiv(e, e2)
				e = div
			case itemCarrot:
				switch e2.ExprType() {
				case CONSTANTF:
					pow := NewPowI(e, int(e2.(*ConstantF).F))
					e = pow
				default:
					pow := NewPowE(e, e2)
					e = pow
				}
			}
			// fmt.Printf("%snew e: %v\n", prefix, e)
		} else if typ == itemIdentifier && itemPrec[itemMul] >= p {
			L.nextToken = next
			e2 := parseExpr(prefix+"  ", L, itemPrec[itemMul])
			mul := NewMul()
			mul.Insert(e)
			mul.Insert(e2)
			e = mul
			// fmt.Printf("%sid-e: %v\ne2: %v\n", prefix, e, e2)
		} else if typ > itemKeyword && itemPrec[itemKeyword] >= p {
			// fmt.Printf("%skey-e: %v\n", prefix, e)
			L.nextToken = next
			e2 := parseExpr(prefix+"  ", L, itemPrec[itemKeyword])
			mul := NewMul()
			mul.Insert(e)
			mul.Insert(e2)
			e = mul
			// fmt.Printf("%se: %v\ne2: %v\n", prefix, e, e2)
		} else if typ == itemLParen || typ == itemLBrack {
			// consumed '(' or '{'
			e2 := parseExpr(prefix+"  ", L, 0)
			// now consume ')','}'
			var next2 item
			if L.nextToken.typ == itemNil {
				next2 = <-L.items
			} else {
				next2, L.nextToken.typ = L.nextToken, itemNil
			}
			typ2 := next2.typ
			if (typ == itemLParen && typ2 != itemRParen) || (typ == itemLBrack && typ2 != itemRBrack) {
				fmt.Printf("error: expected rhs of %v\n", next2)
			}
			// fmt.Printf("%s(): %v  %v\n", prefix, e, e2)

			mul := NewMul()
			mul.Insert(e)
			mul.Insert(e2)
			e = mul
		} else {
			// fmt.Printf("%selse(%d-%d): %v   %v\n", prefix, p, itemPrec[typ], e, next)
			L.nextToken = next
			break
		}
	}
	// fmt.Printf("%sout-E(%d): %v   %v\n", prefix, p, e, L.nextToken)
	return e
}

func parsePiece(prefix string, L *lexer) Expr {
	// fmt.Printf("%sin-P()\n", prefix)
	var next item
	if L.nextToken.typ == itemNil {
		next = <-L.items
	} else {
		next, L.nextToken.typ = L.nextToken, itemNil
	}
	typ := next.typ
	var e Expr
	// fmt.Printf("%spnxt  %v\n", prefix, next)

	// if isUnary(typ) {
	// 	e1 := parseExpr(prefix+"  ", L, 0)
	// 	e = NewNeg(e1)
	// } else 
	if typ == itemNeg {
		e1 := parsePiece(prefix+"  ", L)
		e = NewNeg(e1)
	} else if typ == itemLParen || typ == itemLBrack {
		// consumed '(' or '{'
		e = parseExpr(prefix+"  ", L, 0)
		// now consume ')','}'
		var next2 item
		if L.nextToken.typ == itemNil {
			next2 = <-L.items
		} else {
			next2, L.nextToken.typ = L.nextToken, itemNil
		}
		typ2 := next2.typ
		if (typ == itemLParen && typ2 != itemRParen) || (typ == itemLBrack && typ2 != itemRBrack) {
			fmt.Printf("error: expected rhs of %v\n", next2)
		}
		return e
	} else if typ == itemIdentifier { // leaf
		if next.val[0] == 'C' {
			// coefficient
			ipos := strings.Index(next.val, "_") + 1
			index, err := strconv.ParseInt(next.val[ipos:], 0, 64)
			if err != nil {
				fmt.Printf("Error (%v) parsing index in '%s'\n", err, next.val)
			}
			e = NewConstant(int(index))
		} else if next.val[0] == 'X' {
			// variable
			ipos := strings.Index(next.val, "_") + 1
			index, err := strconv.ParseInt(next.val[ipos:], 0, 64)
			if err != nil {
				fmt.Printf("Error (%v) parsing index in '%s'\n", err, next.val)
			}
			return NewVar(int(index))
		} else {
			// is it a named variable ?
			for p, v := range L.vars {
				if next.val == v {
					e = NewVar(p)
					break
				}
			}
			// is it a named constant ?
			if strings.ToLower(next.val) == "pi" {
				e = NewConstantF(math.Pi)
			}

			if e == nil {
				fmt.Printf("UNIMPLEMENTED IDENTIFIER:  %v\n", next)
			}
		}
	} else if typ == itemNumber {
		flt, err := strconv.ParseFloat(next.val, 64)
		if err != nil {
			fmt.Printf("Error (%v) parsing number in '%s'\n", err, next.val)
		}
		e = NewConstantF(flt)

		// } elsi if { // func name...
	} else if typ > itemKeyword {
		// is it a named function
		switch typ {
		case itemSin:
			e2 := parseExpr(prefix+"  ", L, itemPrec[itemKeyword])
			e = NewSin(e2)
		case itemCos:
			e2 := parseExpr(prefix+"  ", L, itemPrec[itemKeyword])
			e = NewCos(e2)
		case itemTan:
			e2 := parseExpr(prefix+"  ", L, itemPrec[itemKeyword])
			e = NewTan(e2)
		case itemAbs:
			e2 := parseExpr(prefix+"  ", L, itemPrec[itemKeyword])
			e = NewAbs(e2)
		case itemSqrt:
			e2 := parseExpr(prefix+"  ", L, itemPrec[itemKeyword])
			e = NewSqrt(e2)
		case itemExp:
			// now consume the '^'
			var next2 item
			if L.nextToken.typ == itemNil {
				next2 = <-L.items
			} else {
				next2, L.nextToken.typ = L.nextToken, itemNil
			}
			typ2 := next2.typ
			if next.val == "e" && typ2 != itemCarrot {
				fmt.Printf("error: expected ^ after 'e', got %v & %v\n", next, next2)
			}
			e2 := parseExpr(prefix+"  ", L, itemPrec[itemCarrot])
			e = NewExp(e2)
		case itemLog:
			e2 := parseExpr(prefix+"  ", L, itemPrec[itemKeyword])
			e = NewLog(e2)
		}

		if e == nil {
			fmt.Printf("Unimplemented Function:  %v\n", next)
		}

	} else { // error
		L.nextToken = next
		e = nil
	}
	// } else if { is it a user defined function

	// fmt.Printf("%sout-P(): %v\n", prefix, e)
	return e
}
