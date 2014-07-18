#!/bin/bash
# Runs a sample simulation and commits Windows executables, sample output, and docs
# Accepts exactly one argument - version number (like 1.0)
# Should be executed on Windows after finishing all builds

# Test arguments
if [ $# -ne 1 ]; then
  echo "ERROR: requires 1 argument"
  exit 1
fi

ROOTDIR=`pwd`/../
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
svn commit -m "Preparing files for release $1: Windows binaries, sample output, and docs" doc/ sample/ win32/ win64/
if [ $? -ne 0 ]; then
  echo "ERROR: svn commit failed"
  exit 1
fi
