{{define "Title"}}Division{{end}}{{define "PathToStatic"}}../static/{{end}}{{define "Payload"}}
  <header>
    <h1>Division</h1>
  </header>
  <form>
    <p>
      <input type="text" name="dividend" value="{{.Dividend}}" autofocus="true"/> : <input type="text" name="divisor" value="{{.Divisor}}" />
      <input type="submit" value="Submit" />
    </p>
    <p>
      <label for="prec">Precision: </label><input type="text" name="prec" id="prec" value="{{.Precision}}" />
    </p>
    <p>
      <label for="submitprec">Continue calculating until precision, even if remainder is already zero </label><input type="checkbox" id="submitprec" name="stopremz" value="false"{{if not .StopRemz}} checked="checked"{{end}}/>
    </p>
    <p>
      <label for="boxedresult">Display boxes are off </label><input type="checkbox" id="boxedresult" name="boxed" value="false"{{if not .Boxed}} checked="checked"{{end}}/>
    </p>
  </form>{{if .SDivide}}
  <div class="divisionOutputArea">
    {{tplfuncdivdisplay .SDivide .Boxed}}
  </div>{{end}}
  <script type="text/javascript">
    function changeBox() {
      var items = document.getElementsByClassName('divisionColumn');
      for(i=0; i < items.length; i++) {
	if(items[i].hasAttribute('data-division')) {
	  items[i].setAttribute('data-boxed', items[i].getAttribute('data-boxed') == 'true' ? 'false' : 'true');
	}
      }
    }
    document.getElementById('boxedresult').addEventListener('click', changeBox);
    document.getElementById('submitprec').addEventListener('click', function() { document.forms[0].submit();} );
  </script>{{end}}