// Copyright 2011,2012 Johann Höchtl. All rights reserved.
// Use of this source code is governed by a Modified BSD License
// that can be found in the LICENSE file.

// +build !appengine

// RESTFul interface for coordinate transformations - configuration for stand alone server
package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"io/ioutil"
	"os"
)

var configFileName = flag.String("config", "config.json", "location of JSON configuration file")

type config struct {
	APIRoot string
	DocRoot string
	Binding string
	TimeOut int
}

var conf *config

func createorreturnconfig(conf *config) *config {
	if conf == nil {
		conf = &config{APIRoot: "/api", DocRoot: "/doc", Binding: os.Getenv("PORT"), TimeOut: 3600}
	}
	flag.Parse()
	readConfig(*configFileName, conf)
	if conf.Binding == "" {
		conf.Binding = "5000"
	}
	return conf
}

func readConfig(filename string, conf *config) {
	b, err := ioutil.ReadFile(filename)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return
	}

	if conf == nil {
		conf = &config{}
	}

	err = json.Unmarshal(b, &conf)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		panic("Unable to parse json configuration file")
	}
	return
}

func conf_apiroot() string {
	conf = createorreturnconfig(conf)
	return conf.APIRoot
}

func conf_docroot() string {
	conf = createorreturnconfig(conf)
	return conf.DocRoot
}

func conf_binding() string {
	conf = createorreturnconfig(conf)
	return conf.Binding
}

func conf_statictimeout() int {
	conf = createorreturnconfig(conf)
	return conf.TimeOut
}
