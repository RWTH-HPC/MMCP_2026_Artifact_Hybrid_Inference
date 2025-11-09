#!/usr/bin/env python3

CLANG_FORMAT_PATHS = [
    '/bgsys/local/clang/llvm.r201012/bin/clang-format',
    '/pds/opt/llvm/bin/clang-format',
]

import argparse
import itertools
import os
import re
import subprocess
import sys

def main():
    """Main routine that is called if script is directly executed."""

    # Parse command line arguments
    args = parse_arguments()

    # Check if clang-format is on the PATH, otherwise look in known directories
    clang_format_exe = which('clang-format')
    if clang_format_exe is None:
        for fpath in CLANG_FORMAT_PATHS:
            if is_exe(fpath):
                clang_format_exe = fpath
                break
        else:
            abort("Cannot find 'clang-format' executable")


    # Get urls from `svn info`
    svn_url = get_svn_url(args.zfs_root)
    if os.path.basename(svn_url) == 'trunk':
        in_trunk = True
    else:
        in_trunk = False
    svn_root = get_svn_root(args.zfs_root)
    svn_url_trunk = os.path.join(svn_root, 'trunk')


    # Check if this is a branch that has outstanding merges from trunk
    if not in_trunk and has_outstanding_merges(svn_url_trunk, args.zfs_root):
        error("Cannot edit branch that has outstanding merges from trunk.")
        abort("Do a trunk merge (`svn merge ^/trunk`) and try again.")

    # Check if this working copy has local modifications
    if not in_trunk and wc_has_modifications(args.zfs_root) and not args.force:
        error("Cannot edit working copy that has local modifications.")
        if args.build_system:
            abort("Use `commit-force` (instead of `commit`) to override.")
        else:
            abort("Run with '--force' to override.")

    # Obtain diff
    if in_trunk:
        diff = get_svn_diff(os.path.join(args.zfs_root, 'src'))
    else:
        trunk_src_url = os.path.join(svn_url_trunk, 'src')
        diff = get_svn_diff(trunk_src_url, os.path.join(args.zfs_root, 'src'))

    # Parse diff
    hunks = parse_diff(diff)

    # Run clang-format on all found files
    for filename, groups in hunks.items():
        # Skip non-C++ files
        _, ext = os.path.splitext(filename)
        if ext not in ['.h', '.cpp']:
            continue

        # Build & run clang-format command
        command = [clang_format_exe, '-i', '-style=file']
        # Add information on which lines should be formatted
        no_formatted_lines = 0
        for lines in groups:
            command.extend(['-lines', ':'.join([str(l) for l in lines])])
            no_formatted_lines += lines[1] - lines[0] + 1
        # Add filename
        command.append(os.path.join(args.zfs_root, 'src', filename))
        # Run command
        execute(command)

        # Get relative filename
        if os.path.basename(os.getcwd()) == 'src':
            filename_relative = os.path.basename(filename)
        else:
            filename_relative = os.path.join('src', os.path.basename(filename))

        # Print information on formatted lines
        dots = '.' * (80 - 25 - len(filename_relative))
        print("{} {} {:5} line(s) formatted".format(filename_relative, dots,
                no_formatted_lines))

    # Include note to reduce confusion among users about reported changes
    # that do not show up in svn
    print("\nNote: already properly formatted lines are reported but remain "
          "unchanged.")


def has_outstanding_merges(source, target):
    """Check if there are eligible merges from `source` to `target`."""

    out, _, _ = execute(['svn', 'mergeinfo', '--show-revs', 'eligible',
                         source, target])
    return (out != '')


def parse_diff(diff):
    """Parse `diff` as a unified diff and extract files and added lines."""

    hunks = dict()
    lines = []
    filename = None
    in_hunk = False
    line_no = None

    # Parse diff for file and pine information
    for line in diff.splitlines():
        match = re.search('^\+\+\+ (\S+)', line)
        if match:
            if filename is not None:
                if filename != '.':
                    if filename not in hunks:
                        hunks[filename] = []
                    hunks[filename].extend(lines)
                lines = []
            filename = match.group(1)
            in_hunk = False
            continue

        match = re.search('^@@.*\+(\d+),(\d+)', line)
        if match:
            line_no = int(match.group(1))
            in_hunk = True
            continue

        if not line.startswith((' ', '+', '-')):
            in_hunk = False
            continue

        if not in_hunk:
            continue

        if line.startswith('+'):
            lines.append(line_no)
            line_no += 1
        elif line.startswith(' '):
            line_no += 1
        else:
            continue

    # Add last hunk
    if filename != '.':
        if filename not in hunks:
            hunks[filename] = []
        hunks[filename].extend(lines)

    # Delete empty hunks
    hunks = dict((k, v) for k, v in hunks.items() if len(v) > 0)

    # Get line ranges
    for filename, lines in hunks.items():
        groups = []
        for k, g in itertools.groupby(enumerate(lines), lambda x: x[0] - x[1]):
            group = list(map(lambda x: x[1], g))
            groups.append((group[0], group[-1]))
        hunks[filename] = groups

    return hunks


def get_svn_diff(old=None, new=None):
    """Obtain diff from Subversion with `old` and `new` being optional
    arguments.
    """

    command = ['svn', 'diff', '--ignore-properties', '--no-diff-deleted', '--internal-diff']
    if old is not None:
        command.append(old)
        if new is not None:
            command.append(new)
    diff, _, _ = execute(command)
    return diff


def get_svn_url(wc_root):
    """Get the Subversion URL from the given working copy."""

    out, _, _ = execute(['svn', 'info', wc_root])
    for line in out.splitlines():
        if line.startswith('URL:'):
            return line.split()[1]
    return None


def get_svn_root(wc_root):
    """Get the repository root from the given working copy."""

    out, _, _ = execute(['svn', 'info', wc_root])
    for line in out.splitlines():
        if line.startswith('Repository Root:'):
            return line.split()[2]
    return None


def error(msg):
    """Print formatted error message to STDERR  (with colors if writing to a
    terminal).
    """

    formatted = os.path.basename(__file__) + ": error: " + str(msg)
    if sys.stderr.isatty():
        formatted = '\033[31m' + formatted + '\033[0m'
    sys.stderr.write(formatted + '\n')

    # Flush to make sure that formatting codes do not get mixed up
    sys.stderr.flush()


def abort(msg, returncode=1):
    """Print error message and exit with `returncode`."""

    error(msg)
    sys.exit(returncode)


def parse_arguments():
    """Parse and return command line arguments."""

    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('zfs_root', help="Path to ZFS root directory.")
    p.add_argument('--force', action='store_true',
                   help="Continue even if there are local modifications.")
    p.add_argument('--build-system', action='store_true',
                   help="Use flag if invoking from build system.")
    return p.parse_args()


def wc_has_modifications(wc_root):
    """Check if there are uncommitted changes in the given working copy."""

    out, _, _ = execute(['svn', 'st', '-q', wc_root])
    return (out != '')


def execute(command, shell=False, failfast=True):
    """Execute `command`, interpreting it as a shell command if `shell` is true.
    Return (STDOUT, STDERR, returncode) tuple if `failfast` is false or if there
    was no error. If `failfast` is true, STDERR is not piped and the script
    quits immediately upon receiving a non-zero returncode.
    """

    # If failfast is true, do not capture output of stderr
    if failfast:
        stderr = None
    else:
        stderr = subprocess.PIPE

    # Execute command
    p = subprocess.Popen(command, universal_newlines=True, shell=shell,
                         stdout=subprocess.PIPE, stderr=stderr)
    out, err = p.communicate()

    # If failfast is true and the return value indicates an error, abort
    if failfast and p.returncode != 0:
        sys.exit(p.returncode)

    return out, err, p.returncode


def is_exe(path):
    return os.path.isfile(path) and os.access(path, os.X_OK)


def which(executable):
    """Find executable in PATH and return it if found, otherwise return None."""
    for path in os.environ['PATH'].split(':'):
        fpath = os.path.join(path, executable).strip('"')
        if is_exe(fpath):
            return fpath
    else:
        return None


if __name__ == '__main__':
    sys.exit(main())
