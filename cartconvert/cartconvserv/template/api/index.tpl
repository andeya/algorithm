<!DOCTYPE HTML>
<html lang="en">
<head>
  <meta charset="UTF-8"> 
  <title>Cartconvert - API page</title>
  <link rel="shortcut icon" href="/static/images/favicon.png" type="image/png"/> 
  <link rel="icon" href="/static/images/favicon.png" type="image/png"/>
  <link rel="stylesheet" type="text/css" href="/static/css/styles.css"/> 
</head>
<body>
<div>
  <header>
    <h1>Cartconvert - API page</h1>
    <p>Root of API service</p>
  </header>
  <nav>
    <ul>{{range .APIRefs}}
      <li><a href="{{with $.DOCRoot}}{{.}}{{end}}{{.URL}}">{{.Documentation}}</a></li>{{end}}
    </ul>
  </nav>
</div>
</body>
</html>