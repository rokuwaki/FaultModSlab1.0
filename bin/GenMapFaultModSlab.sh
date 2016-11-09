#!/bin/bash 
gmtset MAP_FRAME_TYPE plain
gmtset PS_PAGE_ORIENTATION portrait
PS=knot_map.ps
S=10
line=`gmtinfo knot_value.dat | gawk '{print $8,$7,$9}' | sed -e "s/</ /g" | sed -e "s/>/ /g" | sed -e "s./. .g"`
lon_ll=`echo $line | cut -d ' ' -f 1`
lon_rr=`echo $line | cut -d ' ' -f 2`
lat_bb=`echo $line | cut -d ' ' -f 3`
lat_tt=`echo $line | cut -d ' ' -f 4`
dep_ss=`echo $line | cut -d ' ' -f 5`
dep_dd=`echo $line | cut -d ' ' -f 6`

lon_l=`echo "scale=1; ($lon_ll - 0.5) / 1.0" | bc`
lon_r=`echo "scale=1; ($lon_rr + 0.5) / 1.0" | bc`
lat_b=`echo "scale=1; ($lat_bb - 0.5) / 1.0" | bc`
lat_t=`echo "scale=1; ($lat_tt + 0.5) / 1.0" | bc`
R=$lon_l/$lon_r/$lat_b/$lat_t

dep_d=`echo "scale=0; ($dep_dd + 1) / 1.0" | bc`
dep_s=0
CPT=depth.cpt
makecpt -T$dep_s/$dep_d/1 -Cpolar -Z -I > $CPT

topoCPT=topotopo.cpt
makecpt -T-10000/10000/1000 -Ctopo -Z > $topoCPT

grdcut ~/data/ETOPO1_Bed_g_gmt4.grd -R$R -Gtopo.grd
grdgradient topo.grd -Ggradient.grd -A45
grdimage topo.grd -JM$S -R$R -Igradient.grd -C$topoCPT -Q -t30 -K > $PS
awk '$9=="T" {print $4,$3,$5,$5}' < work/knot_value.dat | psxy -JM$S -R$R -Sc0.05 -C$CPT -K -O >> $PS
psxy ~/data/trench.gmt -J -R -W,gray -O -K >> $PS
pscoast -J -R -W -Ba1f0.5SWen -Df -N1 -O -K >> $PS
psxy work/epicenter.dat -J -R -Sa0.5 -W -: -K -O >> $PS
psscale -C$CPT -D10.5/5/-$S/0.2 -B10/:km: -O >> $PS
psconvert -A -Tf -Z $PS
#ps2raster -A -Tf $PS

rm *.cpt gradient.grd topo.grd
