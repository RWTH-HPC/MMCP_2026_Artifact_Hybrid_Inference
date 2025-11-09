#!/bin/bash

# Load settings
. $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/auxiliary/settings

# Execute regex
$REFAC_TOOL -i spongeBndryId3 `regex_file_list`
