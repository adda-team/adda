#!/bin/bash
# Copyright (C) ADDA contributors
# GNU General Public License version 3
#
# Tests superellipsoid shape in ADDA. Currently relatively simple.
# Should be further developed into a general framework for running equivalent command line combinations in ADDA

# Variables
ADDA='../../src/seq/adda'
DIR1=out1
DIR2=out2
SO1=stdout1
SO2=stdout2
COMMON="-m 1.2 0.01 -size 10 -save_geom"
CLEAN="rm -f -r $SO1 $SO2 $DIR1 $DIR2"
CLEAN2="rm -f -r $SO2 $DIR2"

# Execution function, which adds common options and redirections
function run1 {
  $ADDA $@ $COMMON -dir $DIR1 > $SO1
}
function run2 {
  $ADDA $@ $COMMON -dir $DIR2 > $SO2
}
function run2_compare {
  run2 $@
  diff $DIR1/superellipsoid.geom $DIR2/superellipsoid.geom
  diff $DIR1/CrossSec-X $DIR2/CrossSec-X
  diff $DIR1/CrossSec-Y $DIR2/CrossSec-Y
  $CLEAN2
}

# Sphere
run1 -shape sphere
# Corresponding superellipsoid
run2 -shape superellipsoid 1 1 1 1
# Use diff on resulting geometry files and CrossSec-Y
diff $DIR1/sphere.geom $DIR2/superellipsoid.geom
diff $DIR1/CrossSec-Y $DIR2/CrossSec-Y
$CLEAN

# Ellipsoid
run1 -shape ellipsoid 0.5 2
run2 -shape superellipsoid 0.5 2 1 1
diff $DIR1/ellipsoid.geom $DIR2/superellipsoid.geom
diff $DIR1/CrossSec-X $DIR2/CrossSec-X
diff $DIR1/CrossSec-Y $DIR2/CrossSec-Y
$CLEAN

# Box
run1 -shape box 2 0.5
run2 -shape superellipsoid 2 0.5 0 0
diff $DIR1/box.geom $DIR2/superellipsoid.geom
diff $DIR1/CrossSec-X $DIR2/CrossSec-X
diff $DIR1/CrossSec-Y $DIR2/CrossSec-Y
$CLEAN

# Cylinder
run1 -shape cylinder 1
run2 -shape superellipsoid 1 1 1 0
diff $DIR1/cylinder.geom $DIR2/superellipsoid.geom
diff $DIR1/CrossSec-Y $DIR2/CrossSec-Y
$CLEAN

# Test continuity for small and large arguments
COMMON="-m 1.2 0.01 -size 10 -save_geom -no_vol_cor"
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
COMMON="-m 1.2 0.01 -size 8 -grid 17 -save_geom -no_vol_cor"
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
