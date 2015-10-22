<!DOCTYPE HTML>
<html lang="en">
<head>
  <meta charset="UTF-8"> 
  <title>Cartconvert - Cartography Transformation API</title>
  <link rel="shortcut icon" href="/static/images/favicon.png" type="image/png"/> 
  <link rel="icon" href="/static/images/favicon.png" type="image/png"/>
  <link rel="stylesheet" type="text/css" href="/static/css/styles.css"/> 
</head>
<body>
<div>
  <header>
    <h1>Cartconvert - Cartography Transformation API</h1>
    <p>This service provides a RESTFul API to perform cartography transformations.</p>
  </header>
  <a href="https://github.com/the42/cartconvert">
    <img style="position: absolute; top: 0; right: 0; border: 0;" src="https://s3.amazonaws.com/github/ribbons/forkme_right_orange_ff7600.png" alt="Fork me on GitHub">
  </a>
  <nav>
    <div>
      <ul>
	<li><a href="{{.APIRoot}}">API</a></li>
	{{if .DOCRoot}}<li><a href="{{.DOCRoot.URL}}">{{.DOCRoot.Documentation}}</a></li>{{end}}
      </ul>
    </div>
  </nav>
  <footer>
    <a href="https://twitter.com/myprivate42" class="twitter-follow-button" data-show-count="false" data-size="large">Follow @myprivate42</a>
    <script>!function(d,s,id){var js,fjs=d.getElementsByTagName(s)[0];if(!d.getElementById(id)){js=d.createElement(s);js.id=id;js.src="//platform.twitter.com/widgets.js";fjs.parentNode.insertBefore(js,fjs);}}(document,"script","twitter-wjs");</script>
    <a href="bitcoin:13A8xBk8EEssyNXEba3JzsECmx9fyYVn3D?amount=5&amp;message=Donation%20for%20cartconvert.allowed.org">Like it? Donate!</a>
    <a href="bitcoin:13A8xBk8EEssyNXEba3JzsECmx9fyYVn3D?amount=5&amp;message=Donation%20for%20cartconvert.allowed.org">
      <img class="textmiddle" src="/static/images/accept_bc2s_100x35.png" alt="Bitcoin donate logo"/>
    </a>
  </footer>
</div>
</body>
</html>