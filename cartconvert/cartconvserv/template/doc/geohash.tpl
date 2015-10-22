{{define "Back"}}..{{end}}{{define "Payload"}}
  <header>
    <h1><a href=".">Documentation for Geohash</a></h1>
  </header>
  <h2>Examples</h2>
  <p>
    <a id="osm1" href="#">u4pruydqqvj</a> as <a href="{{.APIRoot}}/geohash/u4pruydqqvj.xml?outputformat=latlongdeg">Lat / Long XML-encoded</a>,
    as <a href="{{.APIRoot}}/geohash/u4pruydqqvj.json?outputformat=utm">UTM, JSON-encoded</a>.
  </p>
  <h2>Reference</h2>
  <p>
    <a href="https://en.wikipedia.org/wiki/Geohash">Wikipedia [EN]</a>, <a href="https://github.com/kungfoo/geohash-java/blob/master/src/main/java/ch/hsr/geohash/GeoHash.java">Java</a>, <a href="http://blog.dixo.net/downloads/geohash-php-class/">PHP</a>
  </p>
  <h2>Geohash API Documentation</h2>
  <p><a href="https://github.com/the42/cartconvert/blob/master/cartconvserv/README.md#geohash---conversions-">Documentation on Github</a> (authorative developer source)
  </p>
  <script>
    document.getElementById("osm1").addEventListener('click', function() {return osmload('{{.APIRoot}}/geohash/u4pruydqqvj.json?outputformat=latlongcomma')});
  </script>
  {{end}}