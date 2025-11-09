#!/bin/bash

# Load settings
. $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/auxiliary/settings

# Execute regex
echo "Qualified properties:"
echo "`regex_file_list`" | xargs -n1 -P10 $SHELL_REGEX_DIR/qualified_properties.sh
echo
echo "Accessors:"
echo "`regex_file_list`" | xargs -n1 -P10 $SHELL_DRIVER b_properties.sh
