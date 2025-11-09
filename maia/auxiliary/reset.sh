#!/bin/bash

# This script 'hopefully' removes all files that were generated during the CMake
# build process, except for the log file of configure.py itself

# Find all build directories (directories containing a Makefile)
# To not mess up legacy code, /src is exluded here
# Documentation is build 'in-source', hence /doc is excluded
build_dirs=$(find . -mindepth 2 -maxdepth 2 -type f -name 'Makefile' -not -path "./doc/*" -not -path "./src/*" -exec dirname '{}' \; | sed 's/\.\///g' | echo $(tr '\012' ' '))

# Set files/directories to remove
to_remove="CMakeCache.txt CMakeFiles Makefile.O cmake_install.cmake \
           Doxyfile build.ninja rules.ninja .ninja_deps .ninja_log
           analyze coverage compile_commands.json ${build_dirs}"

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
