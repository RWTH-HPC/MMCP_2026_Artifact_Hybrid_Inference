#!/bin/bash

# Load settings
. $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/settings

# Check & save variable name
if [ $# -lt 1 ]; then
  echo "error: missing argument: variable name(s)" >&2
  exit 2
fi

# Add underscores to all variables on command line
for v in $*; do
  # Add underscore in cells
  for f in $CELL_FILES; do
    echo "Adding underscore to $v in $SRC_DIR/${f}..."
    grep -q -w $v $SRC_DIR/$f && sed -i "s/\b$v\b/\0_/g" $SRC_DIR/$f
  done

  # Add underscore to accessors
  for f in $ACCESSOR_FILES; do
    echo "Adding underscore to $v in $SRC_DIR/${f}..."
    grep -q "/\*\*/" $SRC_DIR/$f && sed -i "s/\(.*\/\*\*\/.*\)\b$v\b/\1${v}_/g" $SRC_DIR/$f
  done
done

echo "Done."
