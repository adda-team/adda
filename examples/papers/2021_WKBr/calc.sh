#!/bin/bash

# Executes ADDA for a sphere with different approximations for internal fields. Default input parameters lead to 
# lightweight sample simulation. Realistic cases from the paper (Figs.15,16) are given by correponding wrapper scripts,
# but they require a lot more resources (see also the scripts in cluster/).
#
# For full operation, it is required to provide reference and approximate internal fields, using the scripts in WKBr/
# and Mie_solution/ folders. They can be called automatically by this script, but then

# Input parameters. Either change the default values here or specify them in the command line
size=${1:-10}
m=${2:-1.1}
grid=${3:-16}
calcref=${4:-yes}

# Change 4 to the actual number of processor cores to employ
MPIEXEC="mpiexec -n 4"
# Relative path to adda_mpi after default compilation (make mpi)
ADDA_MPI="../../../src/mpi/adda_mpi"
RUN="$MPIEXEC $ADDA_MPI"

# the following define filenames (with full or relative paths) 
EXACT="Mie_solution/if-$size-$m-$grid.dat"
WKBR_I="WKBr/simple_wkbr-$size-$m-$grid-Y-component.dat"
WKBR_II="WKBr/complex_wkbr-$size-$m-$grid-Y-component.dat"

# calculate references (uses Python)
if [ "$calcref" == "yes" ]; then
  # set python environment
  if command -v python3 &> /dev/null; then
    PYTHON=python3
  else
    if command -v python &> /dev/null; then
      PYTHON=python
    else
      echo ERROR: Python environment not found
      exit 1
    fi
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

# The following lines correspond to six lines in each figure part
$RUN -shape sphere -size $size -grid $grid -m $m 0 -no_vol_cor -eps 6 -init_field read $EXACT
$RUN -shape sphere -size $size -grid $grid -m $m 0 -no_vol_cor -eps 6 -init_field inc
$RUN -shape sphere -size $size -grid $grid -m $m 0 -no_vol_cor -eps 6 -init_field zero
$RUN -shape sphere -size $size -grid $grid -m $m 0 -no_vol_cor -eps 6 -init_field wkb
$RUN -shape sphere -size $size -grid $grid -m $m 0 -no_vol_cor -eps 6 -init_field read $WKBR_I
$RUN -shape sphere -size $size -grid $grid -m $m 0 -no_vol_cor -eps 6 -init_field read $WKBR_II
