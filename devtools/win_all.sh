#!/bin/bash
# Master script to builds and commit Windows executables, most functionality is inside child scripts
# Accepts exactly one argument - version number (like 1.4.0)
# Should be executed on Windows
#
# Copyright (C) ADDA contributors
# GNU General Public License version 3

# Test arguments
if [ $# -ne 1 ]; then
  echo "ERROR: requires 1 argument"
  exit 1
fi

SHELL=/bin/bash

$SHELL commit_docs $1
if [ $? -ne 0 ]; then
  echo "ERROR: error during commiting docs"
  exit 1
fi
$SHELL build all ../win64
if [ $? -ne 0 ]; then
  echo "ERROR: error during building ADDA"
  exit 1
fi
$SHELL build_misc ../win64/misc
if [ $? -ne 0 ]; then
  echo "ERROR: error during building misc tools"
  exit 1
fi
$SHELL win_commit.sh $1
if [ $? -ne 0 ]; then
  echo "ERROR: error during commit of binaries"
  exit 1
fi
git tag -a v$1 -m "New release $1"
if [ $? -ne 0 ]; then
  echo "ERROR: creating tag failed"
  exit 1
fi