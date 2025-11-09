#!/usr/bin/env python3

"""
Output information to user that some component was not found.
"""

# Import modules
import sys

# Get arguments
notfound = sys.argv[1:]

# Output message
print("The following components were not found: {c}".format(
    c=','.join(notfound)))
