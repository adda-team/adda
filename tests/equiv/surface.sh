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

$ADDAREF -surf 5 1.3 0.1 -dir dirref >outref
$ADDA -surf 5 1.3 0.1 -dir dir1 >out1
$ADDA -surf 4 1 1 0 1.3 0.1 -dir dir2 >out2
$ADDA -surf 5 1 1.3 0.1 1.3 0.1 -dir dir3 >out3
mydiff dirref/CrossSec-Y dir1/CrossSec-Y 1
mydiff dir1/CrossSec-Y dir2/CrossSec-Y 2
mydiff dir1/CrossSec-Y dir3/CrossSec-Y 3
rm -rf dirref outref dir1 dir2 dir3 out1 out2 out3
