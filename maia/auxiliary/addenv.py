#!/usr/bin/env python3

"""
This script takes the first argument to be the Makefile that should be
modified, the remaining arguments as environment variables that should
be added to said Makefile. If '=' is among the variable name arguments,
all current variables are added, except for `module`, since it is a 
a shell function and not a variable. A backup named `Makefile.O` is
created in the directory where the original `Makefile` resides before
the variables are added.
"""

# Import modules
import os
import shutil
import subprocess
import sys

def main():
    # Save arguments
    makefile_path = sys.argv[1]
    varname_arguments = sys.argv[2:]

    # Check if '=' is among the names, as this indicates that we should use
    # all variables
    if '=' in varname_arguments:
        vars = os.environ
    else:
        vars = varname_arguments

    # Create backup of Makefile
    backupfile_path = os.path.join(os.path.dirname(makefile_path), 'Makefile.O')
    shutil.copy(makefile_path, backupfile_path)

    # Get Makefile contents
    with open(makefile_path, 'r') as f:
        makefile = f.readlines()

    # Add explanatory comments
    makefile.insert(2, "\n")
    makefile.insert(3, "# {fn}: environment variables that need to be preserved\n".format(fn=__file__))
    makefile.insert(4, "# Note: this was requested by the respective host file\n")

    # Add variables but skip 'module' command
    line_no = 5
    for var in [v for v in vars if v in os.environ and not is_shell_function(v)]:
        makefile.insert(line_no, "export {v} = {c}\n".format(v=var, c=os.environ[var].replace('\n',' ')))
        line_no = line_no + 1

    # Save new file contents
    with open(makefile_path, 'w') as f:
        f.write(''.join(makefile))


def is_shell_function(name):
    """Return True if `name` is a shell function."""
    if name in os.environ and os.environ[name].startswith('() {'):
        return True
    else:
        return False


# Run this script only if it is called directly
if __name__ == '__main__':
    main()
