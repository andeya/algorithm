{{define "Back"}}..{{end}}{{define "Payload"}}
  <header>
    <h1><a href=".">Documentation for UTM</a></h1>
  </header>
  <h2>Examples</h2>
  <p>
    <a id="osm1" href="#">17T 630084 4833438</a> as <a href="{{.APIRoot}}/utm/17T%20630084%204833438.xml?outputformat=latlongdeg">Lat / Long in degrees, XML-encoded</a>,
    as <a href="{{.APIRoot}}/utm/17T%20630084%204833438.json?outputformat=geohash">Geohash, JSON-encoded</a>.
  </p>
  <h2>Reference</h2>
  <p>
    <a href="http://en.wikipedia.org/wiki/Universal_Transverse_Mercator_coordinate_system">Wikipedia [EN]</a>, <a href="http://de.wikipedia.org/wiki/UTM-Koordinatensystem">Wikipedia [DE]</a>
  </p>
  <h2>UTM API Documentation</h2>
  <p><a href="https://github.com/the42/cartconvert/blob/master/cartconvserv/README.md#utm---conversions-">Documentation on Github</a> (authorative developer source)
  </p>
  <script>
    document.getElementById("osm1").addEventListener('click', function() {return osmload("{{.APIRoot}}/utm/17T%20630084%204833438.json?outputformat=latlongcomma")});
  </script>
  {{end}}