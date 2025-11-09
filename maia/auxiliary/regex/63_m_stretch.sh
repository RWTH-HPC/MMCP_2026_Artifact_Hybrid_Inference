#!/bin/bash

# Load settings
. $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/auxiliary/settings

# Execute regex
echo "`regex_file_list`" | xargs -n1 -P10 $SHELL_DRIVER m_stretch.sh
