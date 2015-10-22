Copyright 2011,2012 Johann Höchtl. All rights reserved.

Use of this source code is governed by a Modified BSD License
that can be found in the LICENSE file.

This package provides a series of functions to deal with
conversion and transformations of coordinates in the Datum Austria

[http://de.wikipedia.org/wiki/Datum_Austria](http://de.wikipedia.org/wiki/Datum_Austria)

and here specifically of the Bundesmeldenetz, the former federal cartographic datum
of Austria. The Bundesmeldenetz is already widely replaced by UTM coordinates but much legacy
data is still encoded in BMN coordinates. Unlike UTM, the BMN uses the Bessel reference ellipsoid
and uses lat0 at Hierro (canary islands), which makes transformations tedious.
For more information see

DE: [http://www.topsoft.at](http://www.topsoft.at/pstrainer/entwicklung/algorithm/karto/oek/austria_oek.htm#bmn)
EN: [http://www.asprs.org](http://www.asprs.org/resources/grids/03-2004-austria.pdf)

Usage is covered by test cases. For installation and further info navigate to the parent package.
