#!/bin/bash

# Load settings
. $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../settings

# Check arguments
if [ $# -lt 1 ]; then
  echo "error: missing argument: script name" >&2
  exit 2
fi
if [ $# -lt 2 ]; then
  echo "error: missing argument: file name(s)" >&2
  exit 2
fi

# Check & save script name
regex_script=$SHELL_REGEX_DIR/$1
if [ ! -f "$regex_script" ]; then
  echo "error: $regex_script does not exist" >&2
  exit 2
fi

# Run regex script with all supplied filenames
shift
for f in $*; do
  echo "Applying $regex_script to $f..."
  $regex_script $f
done
