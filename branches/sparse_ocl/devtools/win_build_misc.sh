#!/bin/bash
# Builds Windows executables of misc tools
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
MISCDIR="$ROOTDIR/misc"
# We assume that compilation is performed on 64 bit system, thus a separate flag is added only for 32-bit case
if [ $1 -eq "32" ]; then
  XFL="EXTRA_FLAGS=-m32"
fi

# Build misc tools
cd "$ROOTDIR/misc"
# update misc
svn update
if [ $? -ne 0 ]; then
  echo "ERROR: svn update of misc/ failed"
  exit 1
fi
for dir in `ls $MISCDIR`; do
  MAKEFILE="$MISCDIR/$dir/Makefile"
  if [ -f "$MAKEFILE" ]; then
    echo "Processing misc/$dir"
    # This is specific for Windows, ideally should be replaced by something like
    # make install
    dest="$WINDIR/misc/$dir"
    if [ ! -d "$dest" ]; then
        mkdir "$dest"
    fi
    cd "$dest"
    make $XFL -f "$MAKEFILE" srcdir="$MISCDIR/$dir"
    if [ $? -ne 0 ]; then
      echo "ERROR: compilation in 'misc/$dir' failed"
      exit 1
    fi
  fi
done
