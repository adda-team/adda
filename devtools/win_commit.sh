#!/bin/bash
# Runs a sample simulation and commits Windows executables and sample output
# Accepts exactly one argument - version number (like 1.4.0)
# Should be executed on Windows after finishing all builds
#
# Copyright (C) ADDA contributors
# GNU General Public License version 3

# Test arguments
if [ $# -ne 1 ]; then
  echo "ERROR: requires 1 argument"
  exit 1
fi

ROOTDIR="`pwd`/../"
EXECDIR="$ROOTDIR/win64"
RUNDIR="run000_sphere_g16_m1.5"
SAMPLEDIR="$ROOTDIR/sample"

echo "Running sample simulation"
# run and copy sample files
cd $EXECDIR
rm -f -r ExpCount "$RUNDIR"
./adda.exe > stdout
if [ $? -ne 0 ]; then
  echo "ERROR: sample adda run failed"
  exit 1
fi
mv stdout "$SAMPLEDIR"
mv $RUNDIR/* "$SAMPLEDIR/$RUNDIR"
rm -r ExpCount "$RUNDIR"

cd "$ROOTDIR"
git commit -m "Preparing files for release $1: Windows binaries and sample output" sample/ win64/
if [ $? -ne 0 ]; then
  echo "ERROR: git commit failed"
  exit 1
fi
