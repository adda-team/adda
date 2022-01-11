#!/bin/bash

# This is an example.
# The script calculates the exact internal electric field using either Bhfield or Scattnlay. 
# You need to download the required project (Bhfield or Scattnlay) separately.

# Input parameters:
size=10
grid=16
m=1.1
type=scattnlay # program type: [bhfield, scattnlay]
# Path to the folder with the exact-script.sh, bhfield-std.exe and fieldnlay.exe, respectively (paths need to be set):
pathmain=C:/Users/konstantin/Documents/main-script-test/
pathbhfield=C:/Users/konstantin/Documents/bhfield-121005/windows/
pathscattnlay=C:/Users/konstantin/Documents/scattnlay-master/scattnlay-master/
# Calculating the required variables:
rb=$($(which python) find-rborder.py $size $grid)
lb=$($(which python) find-lborder.py $size $grid)
radius=$($(which python) find-radius.py $size)
if [[ "$type" == "bhfield" ]]
then
cd $pathbhfield
bhfield-std.exe 6.28318530718 $radius $radius $grid $lb $rb $grid $lb $rb $grid $lb $rb other 0 1 $m 0 $m 0
cp V_0Ereim.dat $pathmain/bhfield-$size-$m-$grid.dat
elif [[ "$type" == "scattnlay" ]]
then
cd $pathscattnlay
fieldnlay.exe -l 1 $radius $m 0 -p $lb $rb $grid $lb $rb $grid $lb $rb $grid > scattnlay-$size-$m-$grid.dat
cp scattnlay-$size-$m-$grid.dat $pathmain/scattnlay-$size-$m-$grid.dat
fi
