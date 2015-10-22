// Copyright 2012  Johann Höchtl. All rights reserved.
// Use of this source code is governed by a Modified BSD License
// that can be found in the LICENSE file.
//
// configuration for webzapfen
package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"io/ioutil"
	"os"
	"sort"
)

var configFileName = flag.String("config", "config.json", "location of JSON configuration file")

type config struct {
	RootDomain      string
	LocalBinding    string
	OuterPort       string
	RootTemplateDir string
	Languages       map[string]string
	TimeOut         int
}

var conf *config

func readConfig(filename string, conf *config) {
	b, err := ioutil.ReadFile(filename)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return
	}

	err = json.Unmarshal(b, &conf)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		panic("Unable to parse json configuration file")
	}
	return
}

func createorreturnconfig(conf *config) *config {
	if conf == nil {
		conf = &config{RootDomain: "schoolcalc.hoechtl.org",
			LocalBinding:    os.Getenv("PORT"),
			OuterPort:       ":5000",
			Languages:       map[string]string{"de": "Deutsch", "en": "Englisch"},
			RootTemplateDir: "./templates/",
			TimeOut:         3600}

		flag.Parse()
		readConfig(*configFileName, conf)
		if conf.LocalBinding == "" {
			conf.LocalBinding = "5000"
		}
	}
	return conf
}

func conf_rootdomain() string {
	conf = createorreturnconfig(conf)
	return conf.RootDomain
}

// this function returns a sorted list of supported langauge abbreviations
// you should use this slice when to iterate over the full list of languages,
// especially if you need the full list of langauges in a repeatedly consistent
// and sorted manner
func conf_ISOlanguages() []string {
	conf = createorreturnconfig(conf)
	languages := []string{}
	for key, _ := range conf.Languages {
		languages = append(languages, key)
	}
	sort.Strings(languages)
	return languages
}

func conf_languages() map[string]string {
	conf = createorreturnconfig(conf)
	return conf.Languages
}

func conf_localbinding() string {
	conf = createorreturnconfig(conf)
	return conf.LocalBinding
}

func conf_outerport() string {
	conf = createorreturnconfig(conf)
	return conf.OuterPort
}

func conf_roottemplatedir() string {
	conf = createorreturnconfig(conf)
	return conf.RootTemplateDir
}

func conf_statictimeout() int {
	conf = createorreturnconfig(conf)
	return conf.TimeOut
}
