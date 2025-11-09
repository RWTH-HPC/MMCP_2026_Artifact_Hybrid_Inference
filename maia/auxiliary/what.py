#!/usr/bin/env python3

import sys
import textwrap

# Set coloring
color_prefix = "\033[34m"
color_suffix = "\033[0m"

# Save command line arguments
compiler = sys.argv[1]
build_type = sys.argv[2]
host_file = sys.argv[3]
compiler_file = sys.argv[4]
executable = sys.argv[5]
flags = sys.argv[6]
if len(sys.argv) > 7:
    include_dirs = sys.argv[7]
else:
    include_dirs = ""
if len(sys.argv) > 8:
    libraries = sys.argv[8]
else:
    libraries = ""
if len(sys.argv) > 9:
    rpath = sys.argv[9]
else:
    rpath = ""

# Create helper to print formatted output
def prettyprint(field, content, width=80, fieldwidth=16):
    if len(content) == 1:
        content_formatted = textwrap.fill(content[0], width=width-fieldwidth,
                                          subsequent_indent=fieldwidth*' ',
                                          break_long_words=False,
                                          break_on_hyphens=False)
    else:
        content_formatted = ('\n'+fieldwidth*' ').join(content)
    out = "{f:{fw}}{p}{c}{s}".format(f=field+':', fw=fieldwidth,
                                     p=color_prefix, s=color_suffix,
                                     c=content_formatted)
    print(out + '\n')

# Print information
prettyprint("Compiler", [compiler])
prettyprint("Build type",[build_type])
prettyprint("Host file", [host_file])
prettyprint("Compiler file", [compiler_file])
prettyprint("Compiler exec", [executable])
prettyprint("Compiler flags", [flags])
prettyprint("Include dirs", include_dirs.split(' '))
prettyprint("Libraries", libraries.split(' '))
# prettyprint("Runtime path", filter(None, rpath.split(':')))
prettyprint("Runtime path", rpath.split(':'))
