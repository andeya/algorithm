{{define "DocSetup"}}<!DOCTYPE HTML>
<html lang="en">
<head>
  <meta charset="UTF-8"> 
  <title>Documentation for the cartconvserv API{{if .ConcreteHeading}} - {{.ConcreteHeading}}{{end}}</title>
  <link rel="shortcut icon" href="/static/images/favicon.png" type="image/png"/> 
  <link rel="icon" href="/static/images/favicon.png" type="image/png"/>
  <link rel="stylesheet" type="text/css" href="/static/css/styles.css"/> 
</head>
<body>
<div>{{template "Payload" .}}
  <h2>Further Documentation</h2>
  <nav>
    <ul>{{range .Navigation}}
      <li><a href="{{$.DocRoot}}{{.URL}}">{{.Documentation}}</a></li>{{end}}
    </ul>
  </nav>
</div>
<script src="/static/js/cartconvserv.js"></script>
</body>
</html>{{end}}