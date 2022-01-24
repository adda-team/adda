#!bin/sh
# Copyright (C) ADDA contributors
# GNU General Public License version 3
# 
# Since the script is run by 'sh ...' from the Makefile to be compatible with both Unix and Windows (MinGW), we
# define its shebang also as 'sh', which means that it is supposed to be compatible with any posix-compliant shell
# like dash on Ubuntu.
# Extracts git hash of current directory (where script is located) and stores it as a C macro, like
# #define GITHASH "c85b54b"
# in a special file (see variables below). Outputs obtained git hash to stdout. 
# If git hash can't be obtained, the file is emptied (or created empty if doesn't exist).
# File update happens only if it will change the content of the file to avoid redundant rebuilds.

file=githash.h
macro=GITHASH

# if git log is not available, the following should silently produce ""
# git log should produce either text (like "Unversioned directory") or hash string (f09cb11a0d36f94fbe6f78b92ef3f294d0049613)
HASH=`git log --pretty=format:'%h' -n 1 2>/dev/null`
if [ "$HASH" != "" ]; then
  line="#define $macro \"$HASH\""
  if [ -s $file ]; then
    if [ "$(cat $file)" != "$line" ]; then
      echo "$line" > $file
    fi
  else
    echo "$line" > $file
  fi
  echo $HASH
else
  # produces empty file, but only if it is not already empty (and existent)
  if [ \( ! -f $file \) -o \( -s $file \) ]; then
    > $file
  fi
fi
