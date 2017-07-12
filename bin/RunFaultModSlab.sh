#!/usr/local/bin/bash
mkdir -p data work

lat0=46.592  #lat for epicenter
lon0=153.266 #lon for epicenter
latp=48.9    #lat for angular velocity
lonp=-71.7   #lon for angular velocity
w=0.750      #angular velocity (DeMets et al., 2010, GJI, doi:10.1111/j.1365-246X.2009.04491.x)
dl=2.0       #knot interval along strike
dk=2.0       #knot interval along dip
nl=151       #number of knots along `strike` (aparent strike for knot distribution)
nk=151       #number of knots along `dip` (aparent dip for knot distribution)
l0=76        #knot number along strike for hypocenter
k0=76        #knot number along dip for hypocenter
strf=225     #aparent strike angle for knot distribution
Region=kur #alu, mex, cas, izu, ker, kur, phi, ryu, van, sco, sol, sam, sum
           #http://earthquake.usgs.gov/data/slab/
#wget https://earthquake.usgs.gov/data/slab/models/$Region"_slab1.0_clip.xyz"    -P data/
#wget https://earthquake.usgs.gov/data/slab/models/$Region"_slab1.0_strclip.xyz" -P data/
#wget https://earthquake.usgs.gov/data/slab/models/$Region"_slab1.0_dipclip.xyz" -P data/
#wget https://earthquake.usgs.gov/data/slab/models/$Region"_top.in"              -P data/
depthdata='"data/'$Region'_slab1.0_clip.xyz"'
strikedata='"data/'$Region'_slab1.0_strclip.xyz"'
dipdata='"data/'$Region'_slab1.0_dipclip.xyz"'

echo $lat0 $lon0 "0"                    > work/input.dat
echo $latp $lonp $w                    >> work/input.dat
echo "2.0  2.0  151  151  76  76  225" >> work/input.dat
echo $depthdata                        >> work/input.dat
echo $strikedata                       >> work/input.dat
echo $dipdata                          >> work/input.dat

gfortran -o bin/FaultModSlab src/geodesic.f90 src/FaultModSlab.f90
bin/FaultModSlab
rm geodesic.mod
