#!/bin/bash
# Master script to builds and commit Windows executables, all functionality is inside child scripts
# Accepts exactly one argument - version number (like 1.0)
# Should be executed on Windows 64-bit

# Test arguments
if [ $# -ne 1 ]; then
  echo "ERROR: requires 1 argument"
  exit 1
fi

SHELL=/bin/bash

$SHELL win_build_adda.sh 64
if [ $? -ne 0 ]; then
  echo "ERROR: error during building 64-bit ADDA"
  exit 1
fi
$SHELL win_build_adda.sh 32
if [ $? -ne 0 ]; then
  echo "ERROR: error during building 32-bit ADDA"
  exit 1
fi
$SHELL win_build_misc.sh 64
if [ $? -ne 0 ]; then
  echo "ERROR: error during building 64-bit misc tools"
  exit 1
fi
$SHELL win_build_misc.sh 32
if [ $? -ne 0 ]; then
  echo "ERROR: error during building 32-bit misc tools"
  exit 1
fi

$SHELL win_commit.sh $1
if [ $? -ne 0 ]; then
  echo "ERROR: error during commit"
  exit 1
fi
$SHELL make_tag $1
if [ $? -ne 0 ]; then
  echo "ERROR: creating tag failed"
  exit 1
fi