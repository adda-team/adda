#!/bin/bash
# Copyright (C) ADDA contributors
# GNU General Public License version 3
#
# Tests onion (multilayered sphere) and onion_ell (multilayered ellipsoid) shapes in ADDA.

# Variables
ADDA='../../src/seq/adda'
DIR1=out1
DIR2=out2
SO1=stdout1
SO2=stdout2
CLEAN="rm -f -r $SO1 $SO2 $DIR1 $DIR2"
CLEAN2="rm -f -r $SO2 $DIR2"

# Execution function, which adds common options and redirections
function run1 {
  $ADDA $@ $COMMON -dir $DIR1 > $SO1
}
function run2 {
  $ADDA $@ $COMMON -dir $DIR2 > $SO2
}

# sphere and onion (all layers same index)
COMMON="-size 10"
run1 -shape sphere -m 1.2 0.01
run2 -shape onion 0.5 -m 1.2 0.01 1.2 0.01
# Use diff on resulting CrossSec-Y only (geometry files differ)
diff $DIR1/CrossSec-Y $DIR2/CrossSec-Y
$CLEAN

# coated and onion
COMMON="-m 1.2 0.01 1.5 0.1 -size 10 -save_geom"
run1 -shape coated 0.5
run2 -shape onion 0.5
# Use diff on resulting geometry files and CrossSec-Y
diff $DIR1/coated.geom $DIR2/onion.geom
diff $DIR1/CrossSec-Y $DIR2/CrossSec-Y
$CLEAN

# coated2 and onion
COMMON="-m 1.2 0.01 1.5 0.1 1.3 0.001 -size 10 -save_geom"
run1 -shape coated2 0.6 0.2
run2 -shape onion 0.6 0.2
diff $DIR1/coated2.geom $DIR2/onion.geom
diff $DIR1/CrossSec-Y $DIR2/CrossSec-Y
$CLEAN

# ellipsoid and onion_ell (all layers same index)
COMMON="-size 10"
run1 -shape ellipsoid 1.5 0.6 -m 1.2 0.01
run2 -shape onion_ell 1.5 0.6 0.5 -m 1.2 0.01 1.2 0.01
# Use diff on resulting CrossSec-Y only (geometry files differ)
diff $DIR1/CrossSec-Y $DIR2/CrossSec-Y
$CLEAN
