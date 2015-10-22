// +build !nodoc
// Copyright 2011,2012 Johann Höchtl. All rights reserved.
// Use of this source code is governed by a Modified BSD License
// that can be found in the LICENSE file.

// RESTFul interface for coordinate transformations - documentation part
package main

import (
	"bytes"
	"fmt"
	"html/template"
	"log"
	"net/http"
	"net/url"
	"path"
	"runtime/debug"
	"strconv"
)

// These constants specifiy the directory in which the documentation files are saved
const (
	docdefaultsetup = "default.tpl"
	docmainTemplate = "index.tpl" // The main documentation file. Other filenames are created from the requested API documentation
	doctemplateroot = templateroot + "doc/"
)

// defines the layout of a documentation page and is used by html/template
type docPageLayout struct {
	ConcreteHeading  string // Heading for detailed documentation pages
	APIRoot, DocRoot string // APIRoot is used for inline examples
	Navigation       []Link // Keeps an url and the text for a href
}

// user defined APIRoot and DocRoot are constant throughout program execution
var docPage = docPageLayout{APIRoot: conf_apiroot(), DocRoot: conf_docroot()}

func docHandler(w http.ResponseWriter, req *http.Request) {

	// Error handler for documentation
	defer func() {
		if err := recover(); err != nil {
			buf := fmt.Sprintf(httperrorstr, err)
			log.Printf("%s", buf)
			log.Printf("%s", debug.Stack())
			http.Error(w, buf, http.StatusInternalServerError)
		}
	}()

	base := path.Base(req.URL.Path)
	var filename string

	// check if the incoming url is the base url for documentation
	if base == path.Base(docPage.DocRoot) {
		// if the incoming url is the base url for documentation, load the generic help template
		filename = doctemplateroot + docmainTemplate
	} else {
		// else load the specific help template. The filename is constructed from the API function
		filename = doctemplateroot + base + ".tpl"
		docPage.ConcreteHeading = httphandlerfuncs[base].docstring
	}

	tpl, err := template.ParseFiles(filename, doctemplateroot+docdefaultsetup)
	if err != nil {
		panic(err)
	}

	buf := new(bytes.Buffer)
	err = tpl.ExecuteTemplate(buf, "DocSetup", docPage)
	if err != nil {
		panic(err)
	}

	w.Header().Set("Content-Length", strconv.Itoa(buf.Len()))
	buf.WriteTo(w)
}

func init() {
	// parse all REST handlers and create corresponding documentation links
	for _, val := range httphandlerfuncs {
		url, err := url.Parse(val.method + "/")
		if err != nil {
			panic(fmt.Sprintf("%s: %s is not a valid url", err.Error(), val.method))
		}
		docitem := Link{URL: url, Documentation: val.docstring}
		docPage.Navigation = append(docPage.Navigation, docitem)
	}

	// if documentation is compiled in, we want it included as a link on the main page

	url, _ := url.Parse(docPage.DocRoot)
	docrootLink = &Link{URL: url, Documentation: "Documentation"}
	http.HandleFunc(url.String()+"/", docHandler)
}
