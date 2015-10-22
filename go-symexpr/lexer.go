package symexpr

import (
	"fmt"
	"strings"
	"unicode"
	"unicode/utf8"
)

func (l *lexer) run() {
	for state := lexExpr; state != nil; {
		state = state(l)
	}
	l.items <- item{itemEOF, "eof"}
	close(l.items) // No more tokens will be delivered.
}

type item struct {
	typ itemType
	val string
}

func (i item) String() string {
	switch i.typ {
	case itemEOF:
		return "EOF"
	case itemError:
		return i.val
	}
	if len(i.val) > 10 {
		return fmt.Sprintf("%.10q...", i.val)
	}
	return fmt.Sprintf("%d-'%s'", i.typ, i.val)
}

type itemType int

const (
	itemNil itemType = iota
	itemError
	itemEOF // EOF

	itemAdd    // +
	itemMul    // *
	itemDiv    // /
	itemDiv2   // {...}//{...}
	itemNeg    // -
	itemLParen // (
	itemRParen // )
	itemLBrack // {
	itemRBrack // }
	itemCarrot // ^

	itemNumber
	itemIndex
	itemIdentifier // vars,coeffs,funcs

	itemKeyword // placeholder 
	itemFrac    // \frac
	itemCDot    // \cdot  (latex multiplication)
	itemSin
	itemCos
	itemTan
	itemAbs
	itemSqrt
	itemExp
	itemLog
)

func isLeaf(it itemType) bool {
	switch it {
	case itemNumber, itemIdentifier:
		return true
	}
	return false
}
func isUnary(it itemType) bool {
	// switch it {
	// case itemNeg:
	// 	return true
	// }
	return false
}
func isBinary(it itemType) bool {
	switch it {
	case itemAdd, itemNeg, itemMul, itemDiv, itemDiv2, itemCarrot:
		return true
	}
	return false
}

var itemName = map[itemType]string{
	itemNil:   "nil",
	itemError: "error",
	itemEOF:   "EOF",

	itemAdd:    "add",
	itemMul:    "mul",
	itemDiv:    "div",
	itemDiv2:   "div2",
	itemNeg:    "neg",
	itemLParen: "lParen",
	itemRParen: "rParen",
	itemLBrack: "lBrack",
	itemRBrack: "rBrack",
	itemCarrot: "carrot",

	itemNumber:     "number",
	itemIndex:      "index",
	itemIdentifier: "identifier",

	// keywords
	itemCDot: "cdot",
	itemFrac: "frac",
	itemSin:  "sin",
	itemCos:  "cos",
	itemTan:  "tan",
	itemAbs:  "abs",
	itemSqrt: "sqrt",
	itemExp:  "exp",
	itemLog:  "ln",
}

func (i itemType) String() string {
	s := itemName[i]
	if s == "" {
		return fmt.Sprintf("item%d", int(i))
	}
	return s
}

var key = map[string]itemType{
	"\\cdot": itemCDot,
	"\\frac": itemFrac,
	"sin":    itemSin,
	"Sin":    itemSin,
	"cos":    itemCos,
	"Cos":    itemCos,
	"tan":    itemTan,
	"Tan":    itemTan,
	"abs":    itemAbs,
	"Abs":    itemAbs,
	"sqrt":   itemSqrt,
	"Sqrt":   itemSqrt,
	"e":      itemExp,
	"exp":    itemExp,
	"Exp":    itemExp,
	"ln":     itemLog,
}

const eof = -1

// stateFn represents the state of the scanner as a function that returns the next state.
type stateFn func(*lexer) stateFn

// lexer holds the state of the scanner.
type lexer struct {
	name      string    // used only for error reports.
	input     string    // the string being scanned.
	vars      []string  // allowable variables
	state     stateFn   // the next lexing function to enter.
	start     int       // start position of this item.
	pos       int       // current position in the input.
	width     int       // width of last rune read from input.
	items     chan item // channel of scanned items.
	nextToken item
}

// next returns the next rune in the input.
func (l *lexer) next() (r rune) {
	if l.pos >= len(l.input) {
		l.width = 0
		return eof
	}
	r, l.width = utf8.DecodeRuneInString(l.input[l.pos:])
	l.pos += l.width
	return r
}

// peek returns but does not consume the next rune in the input.
func (l *lexer) peek() rune {
	r := l.next()
	l.backup()
	return r
}

// backup steps back one rune. Can only be called once per call of next.
func (l *lexer) backup() {
	l.pos -= l.width
}

// emit passes an item back to the client.
func (l *lexer) emit(t itemType) {
	l.items <- item{t, l.input[l.start:l.pos]}
	l.start = l.pos
}

// ignore skips over the pending input before this point.
func (l *lexer) ignore() {
	l.start = l.pos
}

// accept consumes the next rune if it's from the valid set.
func (l *lexer) accept(valid string) bool {
	if strings.IndexRune(valid, l.next()) >= 0 {
		return true
	}
	l.backup()
	return false
}

// acceptRun consumes a run of runes from the valid set.
func (l *lexer) acceptRun(valid string) {
	for strings.IndexRune(valid, l.next()) >= 0 {
	}
	l.backup()
}

// lineNumber reports which line we're on. Doing it this way
// means we don't have to worry about peek double counting.
func (l *lexer) lineNumber() int {
	return 1 + strings.Count(l.input[:l.pos], "\n")
}

// error returns an error token and terminates the scan by passing
// back a nil pointer that will be the next state, terminating l.nextItem.
func (l *lexer) errorf(format string, args ...interface{}) stateFn {
	l.items <- item{itemError, fmt.Sprintf(format, args...)}
	return nil
}

// nextItem returns the next item from the input.
func (l *lexer) nextItem() item {
	for {
		select {
		case item := <-l.items:
			return item
		default:
			l.state = l.state(l)
		}
	}
	panic("not reached")
}

// lex creates a new scanner for the input string.
func lex(name, input string, vars []string) *lexer {
	l := &lexer{
		name:  name,
		input: input,
		vars:  vars,
		state: lexExpr,
		items: make(chan item, 2), // Two items of buffering is sufficient for all state functions
	}
	return l
}

/*****  State Functions  *****/

// lexExpr is the top level scanner
func lexExpr(l *lexer) stateFn {
	switch r := l.next(); {
	case r == eof || r == '\n':
		return nil
	case isSpace(r):
		l.ignore()
	case r == '+':
		l.emit(itemAdd)
		return lexExpr
	case r == '*':
		l.emit(itemMul)
		return lexExpr
	case r == '^':
		l.emit(itemCarrot)
		return lexExpr
	case r == '/':
		if l.peek() == '/' {
			l.emit(itemDiv2)
		} else {
			l.emit(itemDiv)
		}
		return lexExpr
	case r == '-':
		l.emit(itemNeg)
		return lexExpr
	case r == '(':
		l.emit(itemLParen)
		return lexExpr
	case r == ')':
		l.emit(itemRParen)
		return lexExpr
	case r == '{':
		l.emit(itemLBrack)
		return lexExpr
	case r == '}':
		l.emit(itemRBrack)
		return lexExpr
	case r == '\\':
		// l.backup()
		return lexIdentifier
	case r == '_':
		l.ignore()
		return lexIndex
	case '0' <= r && r <= '9':
		l.backup()
		return lexNumber
	case isAlphaNumeric(r):
		l.backup()
		return lexIdentifier
	default:
		return l.errorf("unrecognized character in action: %#U", r)
	}
	return lexExpr
}

// lexIdentifier scans an alphanumeric or field.
func lexIdentifier(l *lexer) stateFn {
Loop:
	for {
		switch r := l.next(); {
		case isAlphaNumeric(r):
			// absorb.
		case r == '_':
			// absorb and scanIndex
			digits := "0123456789"
			l.acceptRun(digits)
			l.emit(itemIdentifier)
			return lexExpr
		case r == '.' && (l.input[l.start] == '.' || l.input[l.start] == '$'):
			// field chaining; absorb into one token.
		default:
			l.backup()
			word := l.input[l.start:l.pos]
			if !l.atTerminator() {
				return l.errorf("unexpected character %c", r)
			}
			switch {
			case key[word] > itemKeyword:
				l.emit(key[word])
			default:
				l.emit(itemIdentifier)
			}
			break Loop
		}
	}
	return lexExpr
}

// atTerminator reports whether the input is at valid termination character to
// appear after an identifier. Mostly to catch cases like "$x+2" not being
// acceptable without a space, in case we decide one day to implement
// arithmetic.
func (l *lexer) atTerminator() bool {
	r := l.peek()
	if isSpace(r) {
		return true
	}
	switch r {
	case eof, '+', '-', '(', ')', '{', '}', '/', '*', '^':
		return true
	}
	return false
}

// lexNumber scans a number: decimal, octal, hex, float, or imaginary.  This
// isn't a perfect number scanner - for instance it accepts "." and "0x0.2"
// and "089" - but when it's wrong the input is invalid and the parser (via
// strconv) will notice.
func lexNumber(l *lexer) stateFn {
	if !l.scanNumber() {
		return l.errorf("bad number syntax: %q", l.input[l.start:l.pos])
	} else {
		l.emit(itemNumber)
	}
	return lexExpr
}

func lexIndex(l *lexer) stateFn {
	if !l.scanIndex() {
		return l.errorf("bad index syntax: %q", l.input[l.start:l.pos])
	} else {
		l.emit(itemIndex)
	}
	return lexExpr
}

func (l *lexer) scanIndex() bool {
	// Optional leading sign.
	digits := "0123456789"
	l.acceptRun(digits)

	// Next thing mustn't be alphanumeric.
	// if isAlphaNumeric(l.peek()) {
	// 	l.next()
	// 	return false
	// }
	return true
}

func (l *lexer) scanNumber() bool {
	// Optional leading sign.
	l.accept("+-")
	digits := "0123456789"
	l.acceptRun(digits)
	if l.accept(".") {
		l.acceptRun(digits)
	}
	if l.accept("eE") {
		l.accept("+-")
		l.acceptRun("0123456789")
	}
	// Is it imaginary?
	l.accept("i")
	// Next thing mustn't be alphanumeric.
	// if isAlphaNumeric(l.peek()) {
	// 	l.next()
	// 	return false
	// }
	return true
}

// isSpace reports whether r is a space character.
func isSpace(r rune) bool {
	switch r {
	case ' ', '\t', '\n', '\r':
		return true
	}
	return false
}

// isAlphaNumeric reports whether r is an alphabetic or digit
func isAlphaNumeric(r rune) bool {
	return unicode.IsLetter(r) || unicode.IsDigit(r)
}
