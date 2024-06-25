#!/bin/bash
# Copyright (C) ADDA contributors
# GNU General Public License version 3
#
# Common definitions for tests including equivalent command lines. Should be sourced from actual test scripts.
# TODO: there is space for further development into a general framework, but we may prefer Python for that

# Variables
ADDA='../../src/seq/adda'
DIR1=out1
DIR2=out2
SO1=stdout1
SO2=stdout2
CLEAN="rm -f -r $SO1 $SO2 $DIR1 $DIR2"
CLEAN2="rm -f -r $SO2 $DIR2"

status=0
# use of the following function ensures that status equals to the first non-zero exit code
# Moreover, it exits on first error if STOP_ON_EXIT is defined
function test_exit {
  if [[ -n "$STOP_ON_EXIT" ]] && [[ $1 -ne 0 ]]; then
    echo "ERROR: Difference was found (or other error), see output directories for details"
    exit $1
  fi
  if [ $status -eq 0 ]; then
    status=$1
  fi
}
# Execution function, which adds common options and redirections; COMMON can be set before to avoid repetition
function run1 {
  $ADDA $@ $COMMON -dir $DIR1 > $SO1
  test_exit $?
}
function run2 {
  $ADDA $@ $COMMON -dir $DIR2 > $SO2
  test_exit $?
}
# accepts two shapes names, removes the corresponding comment lines from shape files before comparison
function diff_geom {
  sed -i "/'$1'/d" $DIR1/$1.geom
  sed -i "/'$2'/d" $DIR2/$2.geom
  diff $DIR1/$1.geom $DIR2/$2.geom >&2
  test_exit $?
}
# diffs several files from two runs as is
function diff_verbatim {
  for file in $@; do
     diff $DIR1/$file $DIR2/$file >&2
     test_exit $?
  done
}
