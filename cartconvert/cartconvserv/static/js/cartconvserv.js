function createCORSRequest(method, url) {
  var xhr = new XMLHttpRequest();
  if ("withCredentials" in xhr) {

    // Check if the XMLHttpRequest object has a "withCredentials" property.
    // "withCredentials" only exists on XMLHTTPRequest2 objects.
    xhr.open(method, url, true);

  } else if (typeof XDomainRequest != "undefined") {

    // Otherwise, check if XDomainRequest.
    // XDomainRequest only exists in IE, and is IE's way of making CORS requests.
    xhr = new XDomainRequest();
    xhr.open(method, url);

  } else {

    // Otherwise, CORS is not supported by the browser.
    xhr = null;

  }
  return xhr;
}

function osmload(url) {
  xhr = createCORSRequest("GET", url);
  xhr.onload = function() {
    try {
      response = JSON.parse(xhr.responseText);
    }
    catch(err) {
      alert("Parsing of http response as JSON failed"); 
      return;
    }
    latitude = response.Payload.Lat;
    longitude = response.Payload.Long;

    location.href='http://www.openstreetmap.org/index.html?mlat='+ latitude + '&mlon=' + longitude + '&zoom=15';
  }
  xhr.onerror = function() {
    alert(xhr, xhr.status); 
  }
  xhr.send();
  return false;
}