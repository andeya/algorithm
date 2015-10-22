{{define "Back"}}..{{end}}{{define "Payload"}}
  <header>
    <h1><a href=".">Documentation for Bundesmeldenetz</a></h1>
  </header>
  <h2>Examples</h2>
  <p>
    <a id="osm1" href="#">M34 703168 374510</a> as <a href="{{.APIRoot}}/bmn/M34%20703168%20374510.json?outputformat=latlongdeg">Lat / Long JSON-encoded</a>,
    as <a href="{{.APIRoot}}/bmn/M34%20703168%20374510.xml?outputformat=utm">UTM, XML-encoded</a>.
  </p>
  <h2>Reference</h2>
  <p>
    <a href="http://www.asprs.org/resources/grids/03-2004-austria.pdf">asprs.org [EN]</a>, <a href="http://de.wikipedia.org/wiki/Bundesmeldenetz">Wikipedia [DE]</a>
  </p>
  <h2>BMN API Documentation</h2>
  <p><a href="https://github.com/the42/cartconvert/tree/master/cartconvserv/README.md#bmn---conversions-">Documentation on Github</a> (authorative developer source)
  </p>
  <script>
    document.getElementById("osm1").addEventListener('click', function() {return osmload('{{.APIRoot}}/bmn/M34%20703168%20374510.json?outputformat=latlongcomma')});
  </script>
  {{end}}