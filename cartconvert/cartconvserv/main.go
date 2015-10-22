// +build !appengine
// Copyright 2011,2012 Johann Höchtl. All rights reserved.
// Use of this source code is governed by a Modified BSD License
// that can be found in the LICENSE file.

// RESTFul interface for coordinate transformations.
package main

import (
	"bytes"
	"fmt"
	"html/template"
	"log"
	"net/http"
	"strconv"
	"time"
)

const templateroot = "./template/"
const mainTemplate = "index.tpl"

type logRecord struct {
	http.ResponseWriter

	time                time.Time
	ip, method, rawpath string
	responseBytes       int64
	responseStatus      int
	userAgent, referer  string
	proto               string // "HTTP/1.1"
}

// Wrapper arround DefaultServeMux, inspired by
func Log(handler http.Handler) http.Handler {
	return http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		// TODO: add more information
		logr := &logRecord{time: time.Now().UTC(), ip: r.RemoteAddr}
		// TODO: make the logfile format compatible with eg. apache
		log.Printf("%s %s", logr.time, logr.ip)
		handler.ServeHTTP(w, r)
	})
}

func rootHandler(w http.ResponseWriter, req *http.Request) {
	tpl := template.Must(template.ParseFiles(templateroot + mainTemplate))
	rootpage := &apidocpageLayout{APIRoot: apirootLink, DOCRoot: docrootLink}

	buf := new(bytes.Buffer)
	if err := tpl.Execute(buf, rootpage); err != nil {
		http.Error(w, fmt.Sprint(err), http.StatusInternalServerError)
		return
	}
	w.Header().Set("Content-Length", strconv.Itoa(buf.Len()))
	buf.WriteTo(w)
}

// modelled after https://groups.google.com/d/msg/golang-nuts/n-GjwsDlRco/2l8kpJbAHHwJ
// another good read: http://www.mnot.net/cache_docs/
func maxAgeHandler(seconds int, h http.Handler) http.Handler {
	return http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		w.Header().Add("Cache-Control", fmt.Sprintf("max-age=%d, public, must-revalidate, proxy-revalidate", seconds))
		h.ServeHTTP(w, r)
	})
}

func main() {

	http.HandleFunc("/", rootHandler)
	http.Handle("/static/", maxAgeHandler(conf_statictimeout(), http.StripPrefix("/static/", http.FileServer(http.Dir("static")))))
	http.ListenAndServe(":"+conf_binding(), Log(http.DefaultServeMux))
}
