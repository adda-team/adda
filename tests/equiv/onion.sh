#!/bin/bash
# Copyright (C) ADDA contributors
# GNU General Public License version 3
#
# Tests onion (multilayered sphere) and onion_ell (multilayered ellipsoid) against other shapes

# look inside this script for function and variable definitions
source ./common.sh

# uncomment the following for debugging (where exactly the difference appears)
#STOP_ON_EXIT=1

# runs the second case, compares CrossSec-Y with already existing output of run1, and cleans the output2
function run2_comp_csY {
  run2 $@
  diff_verbatim CrossSec-Y
  $CLEAN2
}
# runs the second case, compares both CrossSec with already existing output of run1, and cleans the output2
function run2_comp_csXY {
  run2 $@
  diff_verbatim CrossSec-X CrossSec-Y
  $CLEAN2
}

# in the following we consider standard size, grid, and refractive index, but shorten simulations with -eps
# COMMON_CORE is used for most simulations in this file, COMMON is used for individual blocks
COMMON_CORE="-size 10 -grid 16 -eps 1"

# sphere and onion (all layers same index), geometry files differ
COMMON="$COMMON_CORE"
run1 -shape sphere -m 1.4 0.01
run2 -shape onion 0.5 -m 1.4 0.01 1.4 0.01
diff_verbatim CrossSec-Y
$CLEAN

# coated and onion
COMMON="$COMMON_CORE -m 1.4 0.01 1.2 0.1 -save_geom"
run1 -shape coated 0.5
run2 -shape onion 0.5
diff_geom coated onion
diff_verbatim CrossSec-Y
$CLEAN

# onion, coated 2 and onion_ell - TODO: remove coated when it is removed from ADDA
COMMON="$COMMON_CORE -m 1.4 0.01 1.3 0.1 1.2 0.01 -save_geom"
run1 -shape onion 0.6 0.2
run2 -shape coated2 0.6 0.2
diff_geom onion coated2 
diff_verbatim CrossSec-Y
$CLEAN2
run2 -shape onion_ell 1 1 0.6 0.2
diff_geom onion onion_ell
diff_verbatim CrossSec-Y
$CLEAN2
# the following tests zero-with shells
#COMMON is updated to allow different refractive indices, but we still compare with run1 above
COMMON="$COMMON_CORE -save_geom"
run2_comp_csY -shape onion 1.0 0.6 0.2 -m 2 0 1.4 0.01 1.3 0.1 1.2 0.01
run2_comp_csY -shape onion 0.6 0.6 0.2 -m 1.4 0.01 2 0 1.3 0.1 1.2 0.01 
run2_comp_csY -shape onion 0.6 0.2 0.2 -m 1.4 0.01 1.3 0.1 2 0 1.2 0.01
run2_comp_csY -shape onion 0.6 0.2 0.0 -m 1.4 0.01 1.3 0.1 1.2 0.01 2 0
$CLEAN

# ellipsoid and onion_ell (all layers same index), geometry file differ
COMMON="$COMMON_CORE"
run1 -shape ellipsoid 1.5 0.6 -m 1.4 0.01
run2 -shape onion_ell 1.5 0.6 0.5 -m 1.4 0.01 1.4 0.01
diff_verbatim CrossSec-X CrossSec-Y
$CLEAN

# zero-with shells for onion_ell
COMMON="$COMMON_CORE -save_geom"
run1 -shape onion_ell 1.5 0.6 0.6 0.2 -m 1.4 0.01 1.3 0.1 1.2 0.01 
run2_comp_csXY -shape onion_ell 1.5 0.6 1.0 0.6 0.2 -m 2 0 1.4 0.01 1.3 0.1 1.2 0.01
run2_comp_csXY -shape onion_ell 1.5 0.6 0.6 0.6 0.2 -m 1.4 0.01 2 0 1.3 0.1 1.2 0.01 
run2_comp_csXY -shape onion_ell 1.5 0.6 0.6 0.2 0.2 -m 1.4 0.01 1.3 0.1 2 0 1.2 0.01
run2_comp_csXY -shape onion_ell 1.5 0.6 0.6 0.2 0.0 -m 1.4 0.01 1.3 0.1 1.2 0.01 2 0
$CLEAN

exit $status
