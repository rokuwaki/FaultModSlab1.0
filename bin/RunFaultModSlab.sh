#!/usr/local/bin/bash

mkdir -p data
Region=kur #alu, mex, cas, izu, ker, kur, phi, ryu, van, sco, sol, sam, sum
wget https://earthquake.usgs.gov/data/slab/models/$Region"_slab1.0_clip.xyz"    -P data/
wget https://earthquake.usgs.gov/data/slab/models/$Region"_slab1.0_strclip.xyz" -P data/
wget https://earthquake.usgs.gov/data/slab/models/$Region"_slab1.0_dipclip.xyz" -P data/
wget https://earthquake.usgs.gov/data/slab/models/$Region"_top.in"              -P data/

lat0=46.592
lon0=153.266

depthdata='"data/kur_slab1.0_clip.xyz"'
strikedata='"data/kur_slab1.0_strclip.xyz"'
dipdata='"data/kur_slab1.0_dipclip.xyz"'

mkdir -p work
echo "46.592  153.266  0"               > work/input.dat
echo "48.9  -71.7  0.750"              >> work/input.dat
echo "2.0  2.0  151  151  76  76  225" >> work/input.dat
echo $depthdata                        >> work/input.dat
echo $strikedata                       >> work/input.dat
echo $dipdata                          >> work/input.dat

#lat lon dummy (must be 0) (epicenter)
#lat lon w (Euler pole : NUVEL-1, DeMets et al., 2010, GJI)
#dl dk nl nk l0 k0 strike of entire fault plane

gfortran -o bin/FaultModSlab src/geodesic.f90 src/FaultModSlab.f90
bin/FaultModSlab
rm *.mod
