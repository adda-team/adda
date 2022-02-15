#!/bin/bash

# Executes ADDA for a sphere with different approximations for internal fields. Default input parameters lead to 
# lightweight sample simulation. Realistic cases from the paper (Figs.15,16) are given by correponding wrapper scripts,
# but they require a lot more resources (see also the scripts in cluster/).
#
# For full operation, it is required to provide reference and approximate internal fields, using the scripts in WKBr/
# and Mie_solution/ folders. They are called automatically by this script.
# 
# This script must be called from its parent folder

# Input parameters. Either change the default values here or specify them in the command line
size=${1:-10}
m=${2:-1.1}
grid=${3:-16}
calcref=${4:-yes}

# Change 4 to the actual number of processor cores to employ
MPIEXEC="mpiexec -n 4"

# The following script will locate adda_mpi (often out of the box), but you can override it by setting ADDA_MPI
# explicitly here (uncomment the following line) or somewhere in the environment (the script will return this variable)
# export ADDA_MPI="../../../src/mpi/adda_mpi"
ADDA_MPI=$(../../find_adda mpi)
if [ $? -ne 0 ]; then
  exit 1
fi

RUN="$MPIEXEC $ADDA_MPI"

# the following define filenames (with full or relative paths) 
EXACT="Mie_solution/if-$size-$m-$grid-X-component.dat" # Scattnlay and Bhfield programs return x-polarized el. field
WKBR_I="WKBr/simple_wkbr-$size-$m-$grid-Y-component.dat"
WKBR_II="WKBr/complex_wkbr-$size-$m-$grid-Y-component.dat"

# calculate references (uses Python)
if [ "$calcref" == "yes" ]; then
  # Set python environment. This may cause problems on Windows, if python is installed through WindowsApps
  # For some reason, such Python installation does not work properly from inside the shell script. The following
  # answer describes, how it can be disabled - https://stackoverflow.com/a/57168165/2633728
  if command -v python3 &> /dev/null; then
    PYTHON=python3
  elif command -v python &> /dev/null; then
    PYTHON=python
  else
    echo ERROR: Python environment not found >&2
    exit 1
  fi
  # compute Mie reference
  cd Mie_solution/
  sh exact-script.sh $size $m $grid scattnlay
  $PYTHON toADDA.py -s $size -m $m -g $grid -t scattnlay
  # compute WKBr internal fields
  cd ../WKBr/
  $PYTHON main.py -s $size -m $m -g $grid -t simple_wkbr
  $PYTHON main.py -s $size -m $m -g $grid -t complex_wkbr
  cd ..
fi

# The following lines correspond to six lines in each figure part.
# -orient 90 0 0 option to rotate the coordinate system, because Scattnlay and Bhfield programs return an x-polarized field, 
# while ADDA works with y-polarization by default.
$RUN -shape sphere -size $size -grid $grid -m $m 0 -no_vol_cor -eps 6 -init_field read $EXACT -orient 90 0 0
$RUN -shape sphere -size $size -grid $grid -m $m 0 -no_vol_cor -eps 6 -init_field inc
$RUN -shape sphere -size $size -grid $grid -m $m 0 -no_vol_cor -eps 6 -init_field zero
$RUN -shape sphere -size $size -grid $grid -m $m 0 -no_vol_cor -eps 6 -init_field wkb
$RUN -shape sphere -size $size -grid $grid -m $m 0 -no_vol_cor -eps 6 -init_field read $WKBR_I
$RUN -shape sphere -size $size -grid $grid -m $m 0 -no_vol_cor -eps 6 -init_field read $WKBR_II
