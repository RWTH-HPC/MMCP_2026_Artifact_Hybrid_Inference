#!/bin/bash

# This script 'hopefully' removes all files that were generated during the CMake
# build process, leaving the build directory practically empty

# Set files/directories to remove
to_remove="CMakeCache.txt CMakeFiles Makefile Makefile.O cmake_install.cmake \
           Doxyfile config.log build.ninja rules.ninja .ninja_deps .ninja_log
           analyze coverage"

printf "Making 'distclean'... "

# Remove all files/directory listed above
for i in $to_remove; do
  find . -not -iwholename '*.svn*' -name "$i" -print0 | xargs -0 rm -rf
done

# Remove html and latex directories in doc/ (if existing)
rm -rf doc/code/doxygen/html
rm -rf doc/code/doxygen/latex

# Furthermore, remove empty directories (these could be left-overs from the
# build process as well)
if find . -not -iwholename '*.svn*' -type d -empty | grep -q '.*'; then
  find . -not -iwholename '*.svn*' -type d -empty -print0 | xargs -0 rmdir
fi

echo "done."
