// Copyright 2011,2012 Johann Höchtl. All rights reserved.
// Use of this source code is governed by a Modified BSD License
// that can be found in the LICENSE file.

// RESTFul interface for coordinate transformations.
package main

import (
	"bytes"
	"encoding/json"
	"encoding/xml"
	"fmt"
	"github.com/henrylee2cn/algorithm/cartconvert/cartconvert"
	"github.com/henrylee2cn/algorithm/cartconvert/cartconvert/bmn"
	"github.com/henrylee2cn/algorithm/cartconvert/cartconvert/osgb36"
	"html/template"
	"log"
	"net/http"
	"net/url"
	"path"
	"runtime/debug"
	"strconv"
)

// supported serialization formats
const (
	JSONFormatSpec = ".json"
	XMLFormatSpec  = ".xml"
)

const apitemplateroot = templateroot + "api/"
const apiTemplate = "index.tpl"

// supported representation/transformation formats
const (
	OutputFormatSpec = "outputformat"

	OFlatlongdeg   = "latlongdeg"
	OFlatlongcomma = "latlongcomma"
	OFgeohash      = "geohash"
	OFUTM          = "utm"
	OFBMN          = "bmn"
	OFOSGB         = "osgb"
)

// Interface type for transparent XML / JSON Encoding
type Encoder interface {
	Encode(v interface{}) error
}

// --------------------------------------------------------------------
// Serialization struct definitions
type (
	URLParameter struct {
		Key    string
		Values []string
	}

	GEOConvertRequest struct {
		Method     string
		Value      string
		Parameters []URLParameter
	}

	GEOConvertResponse struct {
		Status            string
		Code              int
		Error             bool
		GEOConvertRequest *GEOConvertRequest // MIND: GEOConvertRequest is named, because XML and JSON serialization behave differently. An unnamed struct element will NOT be serialized by the XML encoder
		Payload           interface{}
	}

	LatLong struct {
		Lat, Long, Fmt string
		LatLongString  string
	}

	GeoHash struct {
		GeoHash string
	}

	UTMCoord struct {
		UTMCoord  *cartconvert.UTMCoord // MIND: UTMCoord is named, because XML and JSON serialization behave differently. An unnamed struct element will NOT be serialized by the XML encoder
		UTMString string
	}

	BMN struct {
		BMNCoord  *bmn.BMNCoord // MIND: BMNCoord is named, because XML and JSON serialization behave differently. An unnamed struct element will NOT be serialized by the XML encoder
		BMNString string
	}

	OSGB36 struct {
		OSGB36Coord  *osgb36.OSGB36Coord // MIND: OSGB36Coord is named, because XML and JSON serialization behave differently. An unnamed struct element will NOT be serialized by the XML encoder
		OSGB36String string
	}
)

// serialize gets called by the respective handler methods to perform the serialization in the requested output representation
func serialize(latlong *cartconvert.PolarCoord, oformat string) (interface{}, error) {
	var serializestruct interface{}
	var err error

	switch oformat {
	case OFlatlongdeg:
		lat, long := cartconvert.LatLongToString(latlong, cartconvert.LLFdms)
		serializestruct = &LatLong{Lat: lat, Long: long, Fmt: cartconvert.LLFdms.String(), LatLongString: latlong.String()}
	case OFlatlongcomma:
		lat, long := cartconvert.LatLongToString(latlong, cartconvert.LLFdeg)
		serializestruct = &LatLong{Lat: lat, Long: long, Fmt: cartconvert.LLFdeg.String(), LatLongString: latlong.String()}
	case OFgeohash:
		serializestruct = &GeoHash{GeoHash: cartconvert.LatLongToGeoHash(latlong)}
	case OFUTM:
		utm := cartconvert.LatLongToUTM(latlong)
		serializestruct = &UTMCoord{UTMCoord: utm, UTMString: utm.String()}
	case OFBMN:
		var bmnval *bmn.BMNCoord
		bmnval, err = bmn.WGS84LatLongToBMN(latlong, bmn.BMNZoneDet)
		if err == nil {
			serializestruct = &BMN{BMNCoord: bmnval, BMNString: bmnval.String()}
		}
	case OFOSGB:
		var osgb36val *osgb36.OSGB36Coord
		osgb36val, err = osgb36.WGS84LatLongToOSGB36(latlong)
		if err == nil {
			serializestruct = &OSGB36{OSGB36Coord: osgb36val, OSGB36String: osgb36val.String()}
		}
	default:
		err = fmt.Errorf("Unsupported output format: '%s'", oformat)
	}
	return serializestruct, err
}

func getfirstValueFromURLParameters(params []URLParameter, key string) (retval string) {
	for _, parameter := range params {
		if parameter.Key == key {
			return parameter.Values[0]
		}
	}
	return
}

// --------------------------------------------------------------------
// http handler methods corresponding to the restful methods
//
func latlongHandler(request *GEOConvertRequest, latlongstrval, oformat string) (interface{}, error) {

	if len(latlongstrval) > 0 {
		return nil, fmt.Errorf("Latlong doesn't accept an input value. Use the parameters 'lat' and 'long' instead")
	}

	slat := getfirstValueFromURLParameters(request.Parameters, "lat")
	slong := getfirstValueFromURLParameters(request.Parameters, "long")

	var lat, long float64
	var err error

	lat, err = cartconvert.ADegMMSSToNum(slat)
	if err != nil {
		lat, err = cartconvert.ADegCommaToNum(slat)
		if err != nil {
			return nil, fmt.Errorf("Not a bearing: '%s'", slat)
		}
	}

	long, err = cartconvert.ADegMMSSToNum(slong)
	if err != nil {
		long, err = cartconvert.ADegCommaToNum(slong)
		if err != nil {
			return nil, fmt.Errorf("Not a bearing: '%s'", slong)
		}
	}

	latlong := &cartconvert.PolarCoord{Latitude: lat, Longitude: long, El: cartconvert.DefaultEllipsoid}
	return serialize(latlong, oformat)
}

func geohashHandler(request *GEOConvertRequest, geohashstrval, oformat string) (interface{}, error) {
	var latlong *cartconvert.PolarCoord
	var err error
	if latlong, err = cartconvert.GeoHashToLatLong(geohashstrval, nil); err != nil {
		return nil, err
	}
	return serialize(latlong, oformat)
}

func utmHandler(req *GEOConvertRequest, utmstrval, oformat string) (interface{}, error) {
	var utmval *cartconvert.UTMCoord
	var err error
	if utmval, err = cartconvert.AUTMToStruct(utmstrval, nil); err != nil {
		return nil, err
	}

	var latlong *cartconvert.PolarCoord
	if latlong, err = cartconvert.UTMToLatLong(utmval); err != nil {
		return nil, err
	}
	return serialize(latlong, oformat)
}

func bmnHandler(req *GEOConvertRequest, bmnstrval, oformat string) (interface{}, error) {
	var bmnval *bmn.BMNCoord
	var err error
	if bmnval, err = bmn.ABMNToStruct(bmnstrval); err != nil {
		return nil, err
	}

	var latlong *cartconvert.PolarCoord
	if latlong, err = bmn.BMNToWGS84LatLong(bmnval); err != nil {
		return nil, err
	}
	return serialize(latlong, oformat)
}

func osgbHandler(req *GEOConvertRequest, osgb36strval, oformat string) (interface{}, error) {
	var osgb36val *osgb36.OSGB36Coord
	var err error
	if osgb36val, err = osgb36.AOSGB36ToStruct(osgb36strval, osgb36.OSGB36Leave); err != nil {
		return nil, err
	}
	return serialize(osgb36.OSGB36ToWGS84LatLong(osgb36val), oformat)
}

// closure of the restful methods
//    enc: requested encoding scheme
//    req: calling context
//    value: coordinate value to be transformed
//    oformat: requested transformation representation, eg. utm, geohash
type restHandler func(resp *GEOConvertRequest, value, oformat string) (interface{}, error)

const httperrorstr = "An error occurred: %s"

func (fn httphandlerfunc) ServeHTTP(w http.ResponseWriter, req *http.Request) {

	// API error handler
	// Recover from panic by setting http error 500 and letting the user know the reason
	defer func() {
		if err := recover(); err != nil {
			buf := fmt.Sprintf(httperrorstr, err)
			log.Printf("%s", buf)
			log.Printf("%s", debug.Stack())
			http.Error(w, buf, http.StatusInternalServerError)
		}
	}()

	if req.ParseForm() != nil {
		panic("Cannot parse request parameters")
	}

	// val: coordinate value
	// serialformat: serialization format
	// oformat: requested output format
	val := path.Base(req.URL.Path)
	serialformat := path.Ext(val)
	val = val[:len(val)-len(serialformat)]
	oformat := req.URL.Query().Get(OutputFormatSpec)

	request := &GEOConvertRequest{Method: fn.method, Value: val}
	for key, value := range req.Form {
		request.Parameters = append(request.Parameters, URLParameter{Key: key, Values: value})
	}

	// enc keeps the requested encoding scheme as requested by content negotiation
	var enc Encoder
	// allocate buffer to which the http stream is written, until it gets responded. By doing so we keep the chance to trap errors and respond them to the caller
	buf := new(bytes.Buffer)

	switch serialformat {
	case JSONFormatSpec, "":
		w.Header().Set("Content-Type", "application/json; charset=utf-8")
		enc = json.NewEncoder(buf)
	case XMLFormatSpec:
		w.Header().Set("Content-Type", "text/xml")
		enc = xml.NewEncoder(buf)
	default:
		panic(fmt.Sprintf("Unsupported serialization format: '%s'", serialformat))
	}

	response := &GEOConvertResponse{GEOConvertRequest: request}

	serial, err := fn.restHandler(request, val, oformat)
	response.Payload = serial
	if err != nil {

		// might as well panic(err) but we add some more info
		// we  serialize the error here in the chosen encoding
		response.Error = true
		response.Status = fmt.Sprint(err)
		w.WriteHeader(http.StatusInternalServerError)
	}

	err = enc.Encode(response)

	if err != nil {
		panic(fmt.Sprintf("Unable to encode response: %s", err))
	}

	// Enable CORS and set the allowed domain to those requested by Origin
	if origin := req.Header.Get("Origin"); origin != "" {
		w.Header().Set("Access-Control-Allow-Origin", origin)
	}

	// prevent chunking by explicitely set the content-length
	w.Header().Set("Content-Length", strconv.Itoa(buf.Len()))
	buf.WriteTo(w)
}

// --------------------------------------------------------------------
// http part of the restful service: Templates, cache-variables, ...
//
type Link struct {
	*url.URL
	Documentation string
}

// docrootLinks gets initialized in docserv.go. It is used here to display a welcome screen
// with reference to documentation (if there is any)
var docrootLink *Link
var apirootLink string
var pageCache []*template.Template

type apidocpageLayout struct {
	APIRoot string
	APIRefs []Link
	DOCRoot *Link
}

func apiHandler(w http.ResponseWriter, req *http.Request) {
	tpl := template.Must(template.ParseFiles(apitemplateroot + apiTemplate))
	apipage := &apidocpageLayout{APIRoot: apirootLink, DOCRoot: docrootLink}
	for _, val := range httphandlerfuncs {
		url, _ := url.Parse(val.method)
		linkitem := Link{URL: url, Documentation: val.docstring}
		apipage.APIRefs = append(apipage.APIRefs, linkitem)
	}

	buf := new(bytes.Buffer)
	if err := tpl.Execute(buf, apipage); err != nil {
		http.Error(w, fmt.Sprint(err), http.StatusInternalServerError)
		return
	}
	w.Header().Set("Content-Length", strconv.Itoa(buf.Len()))
	buf.WriteTo(w)
}

// Definition of restful methods: combine API URI with handler method.
// For every API URI,there may be a corresponding documentation URI
type httphandlerfunc struct {
	method string
	restHandler
	docstring string
}

var httphandlerfuncs = map[string]httphandlerfunc{
	"/latlong": {"/latlong", latlongHandler, "Latitude, Longitude"},
	"/geohash": {"/geohash", geohashHandler, "Geohash"},
	"/utm":     {"/utm", utmHandler, "UTM"},
	"/bmn":     {"/bmn", bmnHandler, "AT:Bundesmeldenetz"},
	"/osgb":    {"/osgb", osgbHandler, "UK:OSGB36"},
}

func init() {
	apirootLink = conf_apiroot()
	http.HandleFunc(apirootLink+"/", apiHandler)

	for _, handle := range httphandlerfuncs {
		http.Handle(apirootLink+handle.method+"/", handle)
	}
}
