#!/bin/bash
# Copyright (C) ADDA contributors
# GNU General Public License version 3
#
# Tests superellipsoid against other shapes. May show differences in last digits

# look inside this script for function and variable definitions
source ./common.sh

# uncomment the following for debugging (where exactly the difference appears)
#STOP_ON_EXIT=1

# runs the second case, compares with already existing output of run1, and cleans the output2
function run2_compare {
  run2 $@
  diff_verbatim superellipsoid.geom CrossSec-X CrossSec-Y
  $CLEAN2
}

# in the following we consider moderate size, grid, and refractive index, but shorten simulations with -eps
# COMMON_CORE is used for most simulations in this file, COMMON is used for individual blocks
COMMON_CORE="-m 1.4 0.01 -size 10 -grid 16 -eps 1 -save_geom"
COMMON="$COMMON_CORE"

# Sphere
run1 -shape sphere
run2 -shape superellipsoid 1 1 1 1
diff_geom sphere superellipsoid
diff_verbatim CrossSec-Y
$CLEAN

# Ellipsoid
# this test incurs round-off errors in r_eff, which may cause noticeable differences for sufficiently large m
run1 -shape ellipsoid 0.5 2
run2 -shape superellipsoid 0.5 2 1 1
diff_geom ellipsoid superellipsoid
diff_verbatim CrossSec-X CrossSec-Y
$CLEAN

# Box
run1 -shape box 2 0.5
run2 -shape superellipsoid 2 0.5 0 0
diff_geom box superellipsoid
diff_verbatim CrossSec-X CrossSec-Y
$CLEAN

# Cylinder
run1 -shape cylinder 1
run2 -shape superellipsoid 1 1 1 0
diff_geom cylinder superellipsoid
diff_verbatim CrossSec-Y
$CLEAN

# Test continuity for small and large arguments
COMMON="$COMMON_CORE -no_vol_cor"
run1 -shape superellipsoid 0.5 2 0 0
run2_compare -shape superellipsoid 0.5 2 1e-10 0
run2_compare -shape superellipsoid 0.5 2 0 1e-10
run2_compare -shape superellipsoid 0.5 2 1e-10 1e-10
run2_compare -shape superellipsoid 0.5 2 1e-10 1e-5
run2_compare -shape superellipsoid 0.5 2 1e-5 1e-10
$CLEAN

run1 -shape superellipsoid 2 0.5 1 0
run2_compare -shape superellipsoid 2 0.5 1 1e-10
run2_compare -shape superellipsoid 2 0.5 1 1e-5
$CLEAN

run1 -shape superellipsoid 2 0.5 0 1
run2_compare -shape superellipsoid 2 0.5 1e-10 1
run2_compare -shape superellipsoid 2 0.5 1e-5 1
$CLEAN

#The grid in the following should be 17x35x9 or 17x9x35 (all odd numbers!)
COMMON="-m 1.4 0.01 -size 10 -grid 17 -eps 1 -save_geom -no_vol_cor"
run1 -shape superellipsoid 2.05 0.5 1 10
run2_compare -shape superellipsoid 2.05 0.5 1 1e5 
run2_compare -shape superellipsoid 2.05 0.5 1 1e10
$CLEAN

run1 -shape superellipsoid 2.05 0.5 10 1
run2_compare -shape superellipsoid 2.05 0.5 1e5 1
run2_compare -shape superellipsoid 2.05 0.5 1e10 1
$CLEAN

run1 -shape superellipsoid 0.5 2.05 10 10
run2_compare -shape superellipsoid 0.5 2.05 1e10 10
run2_compare -shape superellipsoid 0.5 2.05 10 1e10
run2_compare -shape superellipsoid 0.5 2.05 1e10 1e10
$CLEAN

exit $status
