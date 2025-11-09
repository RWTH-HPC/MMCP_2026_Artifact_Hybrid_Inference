#!/bin/sh

# This script does the following:
#
# 1) Use svn revert to recursively revert all changes to the files
# 2) Delete all unversioned files
#
# Author: Michael Schlottke
# Date:   2012-10-26
#
# Original author: Richard Hansen
# Original source: http://stackoverflow.com/a/9144984/1329844
#
# USE WITH CARE!!! This script contains an "rm -rf" statement,
# so if something goes wrong, it will go wrong badly!

# Reset options
REVERT_SVN=0

while getopts ":hr" option; do
  case "$option" in
    h)
      echo "Options for cleansvn:"
      echo
      echo "-h"
      echo "  Show this help."
      echo
      echo "-r"
      echo "  Revert svn before deleting unversioned files."
      exit 0
      ;;
    r)
      REVERT_SVN=1
      ;;
  esac
done

# make sure this script exits with a non-zero return value if the
# current directory is not in a svn working directory
svn info &>/dev/null
if [ $? -ne 0 ]; then
  echo "Error: '`pwd`' is not an SVN working directory." >&2
  exit 1
fi

# First revert svn...
if [ $REVERT_SVN -eq 1 ]; then
  svn revert -R ./
fi

# ...then delete all unversioned files

svn status --no-ignore | grep '^[I?]' | cut -c 9- |
# setting IFS to the empty string ensures that any leading or
# trailing whitespace is not trimmed from the filename
while IFS= read -r f; do
    # tell the user which file is being deleted.  use printf
    # instead of echo because different implementations of echo do
    # different things if the arguments begin with hyphens or
    # contain backslashes; the behavior of printf is consistent
    printf '%s\n' "Deleting ${f}..."
    # if rm -rf can't delete the file, something is wrong so bail
    rm -rf "${f}" || exit 1
done
