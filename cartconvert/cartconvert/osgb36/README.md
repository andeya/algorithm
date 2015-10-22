Copyright 2011, 2012 Johann Höchtl. All rights reserved.
Use of this source code is governed by a Modified BSD License
that can be found in the LICENSE file.

This package provides functions to deal with conversion and transformations of coordinates
of the OSGB36 datum, the UK National Grid.

The conversion between WGS84 and OSGB36 uses a simple helmert transformation,
which in the case of OSGB36 inconsistencies, may result in an accuracy not exceeding +/- 5m, and
for cornercases worse than +/- 15m.

If a higher accuracy is required, a set of helmert parameters must be used or the
procedure described at http://www.ordnancesurvey.co.uk/gps/docs/Geomatics_world.pdf.

For further info see [http://gps.ordnancesurvey.co.uk/etrs89geo_natgrid.asp](http://gps.ordnancesurvey.co.uk/etrs89geo_natgrid.asp)

Usage is covered by test cases. For installation and further info navigate to the parent package.