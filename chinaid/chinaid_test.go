package chinaid

import (
	"testing"
)

func Test18(t *testing.T) {
	// make ChinaID18

	addr := [6]byte{}
	for k, v := range "326221" {
		addr[k] = byte(v)
	}
	birth := [8]byte{}
	for k, v := range "19880102" {
		birth[k] = byte(v)
	}
	seq := [3]byte{}
	for k, v := range "124" {
		seq[k] = byte(v)
	}
	id18, ok := Make18(addr, birth, seq)
	t.Logf("Make18(%v): %s", ok, id18)

	// check ChinaID18
	chinaId, err := Parse(id18.String())
	if err != nil {
		t.Fatalf("%v", err)
	}
	t.Logf("id18 format(%%v): %v, sex: %s\n\n", chinaId, chinaId.SexText())
	t.Logf("id18 format(%%+v): %#v\n\n", chinaId)
}

func Test15(t *testing.T) {
	// make ChinaID15

	addr := [6]byte{}
	for k, v := range "326221" {
		addr[k] = byte(v)
	}
	birth := [6]byte{}
	for k, v := range "880102" {
		birth[k] = byte(v)
	}
	seq := [3]byte{}
	for k, v := range "124" {
		seq[k] = byte(v)
	}
	id15, ok := Make15(addr, birth, seq)
	t.Logf("Make15(%v): %s", ok, id15)

	// check ChinaID15
	chinaId, err := Parse(id15.String())
	if err != nil {
		t.Fatalf("%v", err)
	}
	t.Logf("id15 format(%%v): %v, sex: %s\n\n", chinaId, chinaId.SexText())
	t.Logf("id15 format(%%+v): %#v\n\n", chinaId)
}

func TestValidate(t *testing.T) {
	ok := Validate("32622119880102124X")
	if !ok {
		t.Fail()
	}
}

func TestCheckCode(t *testing.T) {
	checkCode, ok := CheckCode("32622119880102124")
	if !ok || checkCode != "X" {
		t.Fatalf("get '%s', expect 'X'", checkCode)
	}

	checkCode, ok = CheckCode("326221880102124")
	if ok {
		t.Fatalf("get '%s', expect ''", checkCode)
	}
}
