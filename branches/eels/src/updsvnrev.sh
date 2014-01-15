#!bin/sh
# Since the script is run by 'sh ...' from the Makefile to be compatible with both Unix and Windows (MinGW), we
# define its shebang also as 'sh', which means that it is supposed to be compatible with any posix-compliant shell
# like dash on Ubuntu.
# Tests subversion revision number of current directory (where script is located) and stores it as a C macro, like
# #define SVNREV "1234"
# in a special file (see variables below). Outputs obtained revision number to stdout. 
# If revision number can't be obtained, the file is emptied (or created empty if doesn't exist).
# File update happens only if it will change the content of the file to avoid redundant rebuilds.

file=svnrev.h
macro=SVNREV

# The following (commented out) is alternative way to get revision number.
# However, it doesn't test for local modifications), and doesn't seem much faster
#REV=`svn info -R | awk 'BEGIN{max=0} /^Last Changed Rev:/{if ($NF > max) max=$NF} END{printf max}'`

# if svnversion is not available, the following should silently produce ""
# svnversion should produce either text (like "Unversioned directory") or version string (123, 123M, 123:456, 123:456S) 
# awk searches for string starting with number, and removes something like "123:" from the beginning (if present)
REV=`svnversion -c . 2> /dev/null | awk '/^[0-9].*/ {sub(/^[0-9]*:/,"",$1); print $1}'`
if [ "$REV" != "" ]; then
  line="#define $macro \"$REV\""  
  if [ -s $file ]; then
    if [ "$(cat $file)" != "$line" ]; then
      echo "$line" > $file
    fi
  else
    echo "$line" > $file
  fi
  echo $REV
else
  # produces empty file, but only if it is not already empty (and existent)
  if [ \( ! -f $file \) -o \( -s $file \) ]; then
    > $file
  fi
fi
