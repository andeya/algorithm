// Copyright 2011,2012 Johann Höchtl. All rights reserved.
// Use of this source code is governed by a Modified BSD License
// that can be found in the LICENSE file.

// main program for schoolcalc, a web site to help learning to divide.
package main

import (
	"log"
	"net/http"
	"time"
)

type logRecord struct {
	http.ResponseWriter

	time                time.Time
	ip, method, rawpath string
	responseBytes       int64
	responseStatus      int
	userAgent, referer  string
	proto               string // "HTTP/1.1"
}

// Wrapper arround DefaultServeMux
func Log(handler http.Handler) http.Handler {
	return http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		// TODO: add more information
		logr := &logRecord{time: time.Now().UTC(), ip: r.RemoteAddr}
		// TODO: make the logfile format compatible with eg. apache
		log.Printf("%s %s", logr.time, logr.ip)
		handler.ServeHTTP(w, r)
	})
}

func main() {
	http.ListenAndServe(":"+conf_localbinding(), Log(http.DefaultServeMux))
}
