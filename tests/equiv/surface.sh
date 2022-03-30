#!/bin/bash
# Copyright (C) ADDA contributors
# GNU General Public License version 3

# Variables
ADDA='../../src/seq/adda.exe'
ADDAREF='../../win64/adda.exe'

atol=14
rtol=8

function numdiff {
  diff $1 $2 | awk -f ../2exec/diff_numeric.awk -v abs_tol=1e-$atol -v rel_tol=1e-$rtol -  >&2
}

function mydiff {
  tmp=""
  if ! numdiff $1 $2; then
    echo -e "DIFF above is between files '$1' and '$2'" >&2
    echo "TEST $3 FAILED."
  fi

}

$ADDAREF -surf 5 1.3 0.1 -store_beam -dir dirref >outref
$ADDA -surf 5 1.3 0.1 -store_beam -dir dir1 >out1
$ADDA -surf 4 1 0 1 1.3 0.1 -store_beam -dir dir2 >out2
$ADDA -surf 5 1.3 0.1 1 1.3 0.1 -store_beam -dir dir3 >out3
mydiff dirref/IncBeam-Y dir1/IncBeam-Y 1
mydiff dir1/IncBeam-Y dir2/IncBeam-Y 2
mydiff dir1/IncBeam-Y dir3/IncBeam-Y 3
mydiff dirref/CrossSec-Y dir1/CrossSec-Y 4
mydiff dir1/CrossSec-Y dir2/CrossSec-Y 5
mydiff dir1/CrossSec-Y dir3/CrossSec-Y 6
rm -rf dirref outref dir1 dir2 dir3 out1 out2 out3

$ADDAREF -surf 5 1.3 0.1 -store_beam -prop 0 0 -1 -dir dirref >outref
$ADDA -surf 5 1.3 0.1 -store_beam -prop 0 0 -1 -dir dir1 >out1
$ADDA -surf 4 1 0 1 1.3 0.1 -store_beam -prop 0 0 -1 -dir dir2 >out2
$ADDA -surf 5 1.3 0.1 1 1.3 0.1 -store_beam -prop 0 0 -1 -dir dir3 >out3
mydiff dirref/IncBeam-Y dir1/IncBeam-Y 11
mydiff dir1/IncBeam-Y dir2/IncBeam-Y 12
mydiff dir1/IncBeam-Y dir3/IncBeam-Y 13
mydiff dirref/CrossSec-Y dir1/CrossSec-Y 14
mydiff dir1/CrossSec-Y dir2/CrossSec-Y 15
mydiff dir1/CrossSec-Y dir3/CrossSec-Y 16
rm -rf dirref outref dir1 dir2 dir3 out1 out2 out3