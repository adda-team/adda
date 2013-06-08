#!/bin/bash
# Builds Windows executables of main ADDA source
# Accepts exactly one argument - number of bits (32 or 64)
# Should be run on 64-bit Windows

# Test arguments
if [ $# -ne 1 ]; then
  echo "ERROR: requires 1 argument - number of bits for OS (32 or 64)"
  exit 1
fi
if [ $1 -ne "32" -a $1 -ne "64" ]; then
  echo "ERROR: argument 1 should be either 32 or 64"
  exit 1
fi

ROOTDIR=`pwd`/../
WINDIR="$ROOTDIR/win$1"
# We assume that compilation is performed on 64 bit system, thus a separate flag is added only for 32-bit case
if [ $1 -eq "32" ]; then
  XFL="EXTRA_FLAGS=-m32"
fi

# Build ADDA (all versions)
cd "$ROOTDIR/src"
# update src
svn update
if [ $? -ne 0 ]; then
  echo "ERROR: svn update of src/ failed"
  exit 1
fi

# Build FFT versions
make -s $XFL
if [ $? -ne 0 ]; then
  echo "ERROR: compilation of ADDA failed"
  exit 1
fi
# this should be replaced by make install
cp -p seq/adda.exe mpi/adda_mpi.exe ocl/adda_ocl.exe "$WINDIR"

# Build sparse versions
make -s $XFL OPTIONS=SPARSE seq mpi
if [ $? -ne 0 ]; then
  echo "ERROR: compilation of sparse ADDA failed"
  exit 1
fi
# this should be replaced by make install
cp -p seq/adda.exe "$WINDIR"/adda_spa.exe
cp -p mpi/adda_mpi.exe "$WINDIR"/adda_spa_mpi.exe

