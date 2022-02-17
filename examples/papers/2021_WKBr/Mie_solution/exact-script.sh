#!/bin/bash
# The script calculates the exact internal electric field using either Bhfield or Scattnlay. 
# You need to obtain the required executables (Bhfield or Scattnlay) separately and specify the path to the corresponding
# executables below (or copy them in work/ subfolder)

# SET YOUR PARAMETERS (either here or in command line): 
size=${1:-10}
m=${2:-1.1}
grid=${3:-16}
type=${4:-scattnlay} # program type: [bhfield, scattnlay]

# SET YOUR PATHS:
# Path to the working folder (for temporary files, only for bhfield) and for placing output (both must exist)
workpath="$(pwd)/work"
outpath="$(pwd)"
# Executables (with full path) to be used for actual computation
BHFIELD="$workpath/bhfield-std"
# The name of the executable can also be fieldnlay or fieldnlay-mp, depending on installation.
# We tested on the Scattnlay version dated June 7, 2021.
FIELDNLAY="$workpath/fieldnlay-dp"

# Calculating the required variables:
# simple workaround for float computations using awk
calc() { awk "BEGIN{print $*}"; }
radius=$(calc $size/2)
rb=$(calc $radius - $radius/$grid)
lb=-$rb;
if [ "$type" == "bhfield" ]; then
  cd $workpath 
  "$BHFIELD" 6.28318530718 $radius $radius $grid $lb $rb $grid $lb $rb $grid $lb $rb other 0 1 $m 0 $m 0
  cp V_0Ereim.dat "$outpath/bhfield-$size-$m-$grid.dat"
elif [ "$type" == "scattnlay" ]; then
  "$FIELDNLAY" -l 1 $radius $m 0 -p $lb $rb $grid $lb $rb $grid $lb $rb $grid > "$outpath/scattnlay-$size-$m-$grid.dat"
fi
