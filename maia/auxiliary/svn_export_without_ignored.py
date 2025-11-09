#!/usr/bin/env python3

import json
import os
import subprocess
import sys

# Read the canary.ignore JSON file using 'svn cat'
if len(sys.argv) < 4 or len(sys.argv) > 5:
    print('Usage: ' + sys.argv[0] + ' <svnBranch> <canaryIgnore> <testCaseName> [revision] \n'  \
          'svnBranch: URL of the SVN branch to export \n' \
          'canaryIgnore: suffix of the canary.ignore-<canaryIgnore> file \n' \
          'testCaseName: name of the test case \n' \
          'revision: optional SVN revision number to export in the format -r<revision> \n')
    sys.exit(1)
baseUrl = sys.argv[1]
ignoreName = sys.argv[2]
testCaseName = sys.argv[3]
if len(sys.argv) == 5:
    revision = sys.argv[4]
else:
    revision = None
jsonUrl = baseUrl + "/aux/canary.ignore-" + ignoreName + ".json" 
if revision:
    output = subprocess.check_output(['svn', 'cat', revision, jsonUrl]).decode()
else:
    output = subprocess.check_output(['svn', 'cat', jsonUrl]).decode()

# Parse the JSON data and extract the list of testcases to ignore
data = json.loads(output)
ignore_list = data['testcase']['dirs_to_ignore:append']

# Get the base testcase directories with 'svn ls'
if revision:
    output = subprocess.check_output(['svn', 'ls', revision, baseUrl]).decode()
else:
    output = subprocess.check_output(['svn', 'ls', baseUrl]).decode()
lines = output.split('\n')

# Split the directory listing into lines and filter out the ignored directories and files
filtered_lines = []
ignoreListWithSubdirs = [ign.split('/')[0] + '/' for ign in ignore_list if not ign.endswith('/') and '/' in ign]
ignoreListWithSubdirs = list(set(ignoreListWithSubdirs))
for line in lines:
    if not line.strip():
        continue
    if any(ignore_dir in line for ignore_dir in ignore_list):
        print('ignoring ' + line)
        continue
    if line.endswith('/') and line in ignoreListWithSubdirs:
        print('recursing into ' + line)
        suburl = baseUrl + '/' + line
        # Get the subdirectory listing with 'svn ls'
        if revision:
            suboutput = subprocess.check_output(['svn', 'ls', revision, suburl]).decode()
        else:
            suboutput = subprocess.check_output(['svn', 'ls', suburl]).decode()
        sublines = suboutput.split('\n')
        sublines = [line for line in sublines if line.strip()]
        sublines = [line + subline for subline in sublines]
        # Filter out the ignored directories and files
        subfiltered_lines = [line for line in sublines if line.strip() and not any(ignore_dir in line for ignore_dir in ignore_list)]
        filtered_lines.extend([line for line in subfiltered_lines])
    else:
        filtered_lines.append(line)

if not os.path.exists('testcases.' + testCaseName):
    os.makedirs('testcases.' + testCaseName)
for line in filtered_lines:
    if revision:
        print('svn export ' + baseUrl + '/' + line + ' ' + revision + ' testcases.' + testCaseName + '/' + line)
        subprocess.run(['svn' , 'export', baseUrl + '/' + line, revision, 'testcases.' + testCaseName + '/' + line])
    else:
        print('svn export ' + baseUrl + '/' + line + ' testcases.' + testCaseName + '/' + line)
        subprocess.run(['svn' , 'export', baseUrl + '/' + line, 'testcases.' + testCaseName + '/' + line])