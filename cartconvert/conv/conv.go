// Copyright 2011,2012 Johann Höchtl. All rights reserved.
// Use of this source code is governed by a Modified BSD License
// that can be found in the LICENSE file.

// This command reads coordinates from stdin, performs a conversion,
// and writes to stdout. Errors are written to stderr.
//
// The target reference ellipsoid is always the WGS84Ellipsoid
//
// Usage of ./conv
//  -of="deg": specify output format. Possible values are:  dms  geohash  utc  deg
//
package main

import (
	"bufio"
	"flag"
	"fmt"
	"github.com/henrylee2cn/algorithm/cartconvert/cartconvert"
	"github.com/henrylee2cn/algorithm/cartconvert/cartconvert/bmn"
	"github.com/henrylee2cn/algorithm/cartconvert/cartconvert/osgb36"
	"io"
	"os"
	"strings"
)

type displayformat byte

const (
	offmtunknown displayformat = iota
	ofdeg
	ofdms
	ofutm
	ofgeohash
)

type inputformat byte

const (
	ifbmn = iota
	ifosgb36
)

var ofOptions = map[string]displayformat{"deg": ofdeg, "dms": ofdms, "utm": ofutm, "geohash": ofgeohash}
var ifOptions = map[string]inputformat{"bmn": ifbmn, "osgb36": ifosgb36}

func main() {

	var ofcmdlinespec, ifcmdlinespec string
	var of displayformat
	var ifm inputformat
	var lines uint
	var instring, outstring, ofparamvalues, ifparamvalues string
	var pc *cartconvert.PolarCoord

	for key, _ := range ofOptions {
		ofparamvalues += fmt.Sprintf(" %s ", key)
	}

	for key, _ := range ifOptions {
		ifparamvalues += fmt.Sprintf(" %s ", key)
	}

	flag.StringVar(&ofcmdlinespec, "of", "deg", "specify output format. Possible values are: "+ofparamvalues)
	flag.StringVar(&ifcmdlinespec, "if", "osgb36", "specify input format. Possible values are: "+ifparamvalues)
	flag.Parse()

	of = ofOptions[strings.ToLower(ofcmdlinespec)]
	ifm = ifOptions[strings.ToLower(ifcmdlinespec)]

	reader := bufio.NewReaderSize(os.Stdin, 100)
	longline := false

	for data, prefix, err := reader.ReadLine(); err != io.EOF; data, prefix, err = reader.ReadLine() {
		if err != nil {
			fmt.Fprintf(os.Stderr, "conv %d: %s\n", lines, err)
			continue
		}

		if prefix {
			longline = true
			continue
		}

		if longline {
			longline = false
			continue
		}

		lines++

		instring = strings.TrimSpace(string(data))

		if len(instring) == 0 {
			continue
		}

		switch ifm {
		case ifbmn:

			bmncoord, err := bmn.ABMNToStruct(instring)

			if err != nil {
				fmt.Fprintf(os.Stderr, "BMN: error on line %d: %s\n", lines, err)
				continue
			}
			pc, err = bmn.BMNToWGS84LatLong(bmncoord)

			if err != nil {
				fmt.Fprintf(os.Stderr, "BMN: error on line %d: %s (BMN does not return a lat/long bearing)\n", lines, err)
				continue
			}
		case ifosgb36:
			osgb36coord, err := osgb36.AOSGB36ToStruct(instring, osgb36.OSGB36Auto)

			if err != nil {
				fmt.Fprintf(os.Stderr, "OSGB36: error on line %d: %s\n", lines, err)
				continue
			}
			pc = osgb36.OSGB36ToWGS84LatLong(osgb36coord)
		}

		switch of {
		case ofdeg:
			lat, long := cartconvert.LatLongToString(pc, cartconvert.LLFdeg)
			outstring = lat + ", " + long
		case ofdms:
			outstring = pc.String()
		case ofutm:
			outstring = cartconvert.LatLongToUTM(pc).String()
		case ofgeohash:
			outstring = cartconvert.LatLongToGeoHash(pc)
		default:
			fmt.Fprintln(os.Stderr, "Unrecognized output specifier")
			flag.Usage()
			fmt.Fprintf(os.Stderr, "possible values are: [%s]\n", ofparamvalues)
			fmt.Fprintln(os.Stderr, "]")
			os.Exit(2)
		}
		fmt.Fprintf(os.Stdout, "%s\n", outstring)
	}
}
