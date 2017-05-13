package chinaid

import (
	"encoding/json"
	"errors"
	"fmt"
	"strconv"
	"unicode"
	"unsafe"
)

// ErrInvalidChinaID error
var ErrInvalidChinaID = errors.New("invalid Chinese identity card number")

// ChinaID Chinese identity card number
type ChinaID struct {
	whole []byte
	body  []byte
	addr  []byte
	birth []byte
	seq   []byte
	sex   int // 1 male / 0 female
	check []byte
}

// Make18 generates an 18-character length number
func Make18(addr [6]byte, birth [8]byte, seq [3]byte) (*ChinaID, bool) {
	body := string(addr[:]) + string(birth[:]) + string(seq[:])
	if !isDigit(body, -1) {
		return nil, false
	}
	c := &ChinaID{
		body:  []byte(body),
		addr:  addr[:],
		birth: birth[:],
		seq:   seq[:],
	}
	c.check = []byte{c.createCheckCode()}
	c.whole = append(c.body, c.check...)
	c.setSex()
	return c, true
}

// Make15 generates an 15-character length number
func Make15(addr [6]byte, birth [6]byte, seq [3]byte) (*ChinaID, bool) {
	body := string(addr[:]) + string(birth[:]) + string(seq[:])
	if !isDigit(body, -1) {
		return nil, false
	}
	c := &ChinaID{
		body:  []byte(body),
		addr:  addr[:],
		birth: birth[:],
		seq:   seq[:],
	}
	c.whole = c.body
	c.setSex()
	return c, true
}

// Validate judge whether chinese identity card number is valid.
func Validate(id string) bool {
	l := len(id)
	switch l {
	case 15:
		return true
	case 18:
		var sum int
		for k, v := range id[:l-1] {
			d, _ := strconv.Atoi(string(v))
			sum += weight[k] * d
		}
		return validate[sum%11] == id[l-1]
	default:
		return false
	}
}

// CheckCode generates check code of chinese identity card number.
func CheckCode(idBody string) (string, bool) {
	switch len(idBody) {
	case 15:
		return "", true
	case 17:
		var sum int
		for k, v := range idBody {
			d, _ := strconv.Atoi(string(v))
			sum += weight[k] * d
		}
		return string(validate[sum%11]), true
	default:
		return "", false
	}
}

// Parse parse chinese identity card number to *ChinaID
func Parse(id string) (*ChinaID, error) {
	if !isDigit(id, 17) {
		return nil, ErrInvalidChinaID
	}
	whole := []byte(id)
	switch len(whole) {
	case 18:
		c := &ChinaID{
			whole: whole,
			body:  whole[:17],
			addr:  whole[:6],
			birth: whole[6:14],
			seq:   whole[14:17],
			check: whole[17:],
		}
		if !c.validate() {
			return nil, ErrInvalidChinaID
		}
		c.setSex()
		return c, nil

	case 15:
		c := &ChinaID{
			whole: whole,
			body:  whole,
			addr:  whole[:6],
			birth: whole[6:12],
			seq:   whole[12:],
		}
		c.setSex()
		return c, nil
	default:
		return nil, ErrInvalidChinaID
	}
}

// String returns chinese identity card number
func (c *ChinaID) String() string {
	return bytes2String(c.whole)
}

// MarshalJSON implements json.Marshaler
func (c *ChinaID) MarshalJSON() ([]byte, error) {
	var m = map[string]interface{}{
		"whole": c.String(),
		"body":  c.Body(),
		"addr":  c.Addr(),
		"birth": c.Birth(),
		"seq":   c.Seq(),
		"sex":   c.sex,
	}
	if len(c.check) > 0 {
		m["check"] = c.Check()
	}
	return json.Marshal(m)
}

// func (c *ChinaID) GoString() string {
// 	b, _ := c.MarshalJSON()
// 	return bytes2String(b)
// }

// Format implements fmt.Formatter
func (c *ChinaID) Format(f fmt.State, r rune) {
	if r == 'v' && (f.Flag('+') || f.Flag('#')) {
		b, _ := c.MarshalJSON()
		f.Write(b)
		return
	}
	f.Write(c.whole)
}

// Body returns main part of chinese identity card number
func (c *ChinaID) Body() string {
	return bytes2String(c.body)
}

// Addr returns address code
func (c *ChinaID) Addr() string {
	return bytes2String(c.addr)
}

// Birth returns date of birth
func (c *ChinaID) Birth() string {
	return bytes2String(c.birth)
}

// Seq returns sequence number in the persons whose date of birth are same
func (c *ChinaID) Seq() string {
	return bytes2String(c.seq)
}

// Check returns check code
func (c *ChinaID) Check() string {
	return bytes2String(c.check)
}

// Sex returns sex
// 1 male / 0 female
func (c *ChinaID) Sex() int {
	return c.sex
}

// SexText returns sex text, male or female
func (c *ChinaID) SexText() string {
	if c.IsMale() {
		return "male"
	}
	return "female"
}

// IsMale judge whether is male or not
func (c *ChinaID) IsMale() bool {
	return c.sex == 1
}

// IsMale judge whether is female or not
func (c *ChinaID) IsFemale() bool {
	return c.sex == 0
}

func (c *ChinaID) setSex() {
	i, _ := strconv.Atoi(bytes2String(c.seq[2:]))
	c.sex = i % 2
}

var (
	weight   = [17]int{7, 9, 10, 5, 8, 4, 2, 1, 6, 3, 7, 9, 10, 5, 8, 4, 2}
	validate = [17]byte{'1', '0', 'X', '9', '8', '7', '6', '5', '4', '3', '2'}
)

func (c *ChinaID) createCheckCode() byte {
	var sum int
	for k, v := range c.body {
		d, _ := strconv.Atoi(string(v))
		sum += weight[k] * d
	}
	return validate[sum%11]
}

func (c *ChinaID) validate() bool {
	return c.check[0] == c.createCheckCode()
}

func isDigit(s string, maxLen int) bool {
	if maxLen < 0 {
		for _, r := range s {
			if !unicode.IsDigit(r) {
				return false
			}
		}
	} else {
		for i, r := range s {
			if i >= maxLen {
				break
			}
			if !unicode.IsDigit(r) {
				return false
			}
		}
	}
	return true
}

func bytes2String(b []byte) string {
	return *(*string)(unsafe.Pointer(&b))
}
