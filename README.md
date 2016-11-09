# FaultModSlab1.0

- Codes for generating non-planar fault model by interpolating slab1.0 model ([Hayes et al., 2011, JGR](https://earthquake.usgs.gov/data/slab/))
- Conversion between local cartesian and geographic coordinates is performed with [GeographicLib](http://geographiclib.sourceforge.net). Converted (to Fortran90) version of geographiclib is used (provided by [alex-robinson/coordinates](https://github.com/alex-robinson/coordinates))
- Rake angle for each source knot is assumed to be a pure thrust against plate motion derived from MORVEL model ([DeMets et al., 2010, GJI](http://doi.org/10.1111/j.1365-246X.2009.04491.x)).
- The codes are achieved with significant contribution from Amato Kasahara.

## Usage

```bash
$ bash bin/RunFaultModSlab.sh
```

- Before executing script, edit input parameters controlling fault model and slab1.0 models to be considered.

## Dependencies

- GNU bash, version 4.4.0(1)
- GNU Fortran (Homebrew gcc 6.2.0) 6.2.0
- GeographicLib 1.46
