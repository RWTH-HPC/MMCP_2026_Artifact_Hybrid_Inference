#!/bin/bash

# Load settings
. $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/auxiliary/settings

# Run all scripts
all_scripts=$(ls -1 $REGEX_DIR/[0-9][0-9]* | sort)
for s in $all_scripts; do
  echo "Running $s..."
  $s > /dev/null
done
echo "Done."
