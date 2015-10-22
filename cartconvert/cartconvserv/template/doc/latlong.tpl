{{define "Back"}}..{{end}}{{define "Payload"}}
  <header>
    <h1><a href=".">Documentation for Latitude / Longitude</a></h1>
  </header>
  <h2>Examples</h2>
  <p>
    <a id="osm1" href="#">Lat 47.57° Lon 14°0'27''</a> as <a href="{{.APIRoot}}/latlong/.json?lat=47.57°&amp;long=14°0'27''&amp;outputformat=utm">UTM JSON-encoded</a>,
    as <a href="{{.APIRoot}}/latlong/.xml?lat=47.57°&amp;long=14°0'27''&amp;outputformat=latlongcomma">Lat / Long in fractions, XML-encoded</a>.
  </p>
  <h2>Reference</h2>
  <p>
    <a href="http://en.wikipedia.org/wiki/Geographic_coordinate_system">Wikipedia [EN]</a>, <a href="http://en.wikipedia.org/wiki/Geographic_coordinate_system">Wikipedia [DE]</a>
  </p>
  <h2>Lat / Long API Documentation</h2>
  <p><a href="https://github.com/the42/cartconvert/blob/master/cartconvserv/README.md#latitude--longitude---conversions-">Documentation on Github</a> (authorative developer source)
  </p>
  <script>
    document.getElementById("osm1").addEventListener('click', function() {return osmload("{{.APIRoot}}/latlong/.json?lat=47.57°&amp;long=14°0'27''&amp;outputformat=latlongcomma")});
  </script>
  {{end}}