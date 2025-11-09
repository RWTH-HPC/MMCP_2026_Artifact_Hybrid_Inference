#!/usr/bin/env python3

gcov_data_file_extensions = ['.gcno', '.gcda']
gcov_report_file_extensions = ['.gcov']
report_categories = [
    ('DG', ['zfsdg']),
    ('FV', ['zfsfv']),
    ('GEOMETRY', ['zfsgeometry', 'zfsellipsoid']),
    ('GRID', ['zfscartesian', 'zfsgridgen', 'zfspointbased']),
    ('IO', ['zfsinfoout', 'zfsbigdata', 'zfscontext', 'zfsio', 'zfsproperty']),
    ('LBM', ['zfslbm']),
    ('LS', ['zfsls']),
    ('POSTPROCESSING', ['zfsavg', 'zfspostprocessing']),
    ('STRUCTURED', ['zfsstrctrd']),
    ('OTHER', ['']),
]
zfs_file_filter = ['zfs']

import argparse
import os
import shutil
import subprocess
import sys


def main():
    # Get command line arguments
    args = parse_arguments()

    # Determine file name of summary file
    summary_file = os.path.join(args.output_directory, args.summary_file)

    if has_new_coverage_data(args.object_dir, summary_file) or args.force:
        # If new data is available, re-generate the report summary
        generate_reports(args.source_dir, args.object_dir, args.gcov)
        move_reports(args.output_directory)
        data = parse_reports(args.output_directory)
        write_summary_file(summary_file, data)
    else:
        # Read summary file and output to user
        data = read_summary_file(summary_file)

    # Print coverage data
    print_coverage_data(data, args.v, args.output_directory)


def colorize(message, kind):
    if kind == 'GOOD':
        return '\033[32m' + str(message) + '\033[0m'
    elif kind == 'OK':
        return '\033[34m' + str(message) + '\033[0m'
    elif kind == 'BAD':
        return '\033[31m' + str(message) + '\033[0m'
    else:
        raise ValueError("'kind' must be 'GOOD', 'OK', or 'BAD'")


def get_color_kind(percentage):
    if percentage <= 50.0:
        return 'BAD'
    elif percentage <= 90.0:
        return 'OK'
    else:
        return 'GOOD'


def print_coverage_data(data, verbose, output_directory):
    print('')
    print(80 * '=')
    print("ZFS COVERAGE DATA")
    print(80 * '=')

    lines_total = 0
    lines_executed = 0
    for cat, filters in report_categories:
        category_data = []
        total = 0
        executed = 0
        filtered = filter(lambda x: x[0].startswith(tuple(filters)), data)
        for filename, tot, ex in filtered:
            total += tot
            lines_total += tot
            executed += ex
            lines_executed += ex
        data = [x for x in data if x not in filtered]

        if total == 0:
            percentage = 0.0
        else:
            percentage = float(executed)/total * 100
        per = colorize("{:5.1f}%".format(percentage),
                       get_color_kind(percentage))
        dots = (80 - len(cat) - 6 - 2) * '.'
        print("{cat} {dots} {per}".format(cat=cat, dots=dots, per=per))
        if verbose:
            for filename, total, executed in filtered:
                if total == 0:
                    percentage = 0.0
                else:
                    percentage = float(executed)/total * 100
                dots = (80 - 2 - len(filename) - 6 - 2) * '.'
                per = colorize("{:5.1f}%".format(percentage),
                               get_color_kind(percentage))
                print("  {fn} {dots} {per} ({ex:5}/{tot:5} lines)"
                      .format(fn=filename, dots=dots, per=per, ex=executed, tot=total))

            print('')
    print(80 * '-')
    cat = 'TOTAL'
    if lines_total == 0:
        percentage = 0.0
    else:
        percentage = float(lines_executed)/lines_total * 100
    per = colorize("{:5.1f}%".format(percentage), get_color_kind(percentage))
    print("{cat:74}{per}".format(cat=cat, per=per))
    print(80 * '=')
    print('')
    good = colorize("More than 90% is good,", 'GOOD')
    ok = colorize("more than 50% is OK,", 'OK')
    bad = colorize("less than 50% is bad.", 'BAD')
    print("Note: {good} {ok} {bad}".format(bad=bad, good=good, ok=ok))
    print("Note: '{}/' contains detailed information for each file."
          .format(output_directory))
    if verbose:
        print("Note: run 'make/ninja coverage' to get a less detailed report.")
    else:
        print("Note: run 'make/ninja coverage-verbose' "
              "to get a more detailed report.")


def read_summary_file(summary_file):
    data = []
    with open(summary_file, 'rb') as f:
        for line in f:
            filename, lines_total, lines_executed = line.split(':')
            data.append((filename, int(lines_total), int(lines_executed)))
    return data


def write_summary_file(summary_file, data):
    with open(summary_file, 'wb') as f:
        for filename, lines_total, lines_executed in data:
            f.write("{}:{}:{}\n".format(filename, lines_total, lines_executed))


def has_new_coverage_data(object_dir, summary_file):
    if not os.path.isfile(summary_file):
        return True
    last_mod_time = 0
    for f in os.listdir(object_dir):
        base, ext = os.path.splitext(f)
        if ext not in ['.gcda', '.gcno']:
            continue
        path = os.path.join(object_dir, f)
        last_mod_time = max(os.path.getmtime(path), last_mod_time)
    return os.path.getmtime(summary_file) < last_mod_time


def parse_reports(dirname):
    data = []
    for f in sorted(os.listdir(dirname)):
        base, ext = os.path.splitext(f)
        if ext not in gcov_report_file_extensions:
            continue
        path = os.path.join(dirname, f)
        report = parse_gcov_report(path)
        total = report[0]
        executed = len(report[1])
        data.append((base, total, executed))
    return data


def parse_gcov_report(filename):
    lines_executed = []
    lines_total = 0
    with open(filename, 'rb') as f:
        for l in f:
            #skip lines which indicate template specialization
            if '-:' not in l and '#:' not in l:
              tmp = l.split(':', 1)[0].strip()
              number = True
              try:
                tmp = int(tmp)
              except ValueError:
                number = False
              if not number:
                continue

            count, line_no, line = [x.strip() for x in l.split(':', 2)]

            if count == '-':
                continue
            lines_total += 1
            if count not in ['#####', '====']:
                lines_executed.append(line_no)

    return lines_total, lines_executed


def move_reports(dirname, path='.'):
    if not os.path.isdir(dirname):
        os.makedirs(dirname)
    for f in os.listdir(path):
        if f.endswith('.gcov'):
            src = os.path.join(path, f)
            dst = os.path.join(dirname, f)
            if not f.startswith(tuple(zfs_file_filter)):
                os.remove(src)
            else:
                shutil.move(src, dst)


def generate_reports(source_dir, object_dir, gcov_exe='/pds/opt/gcc/bin/gcov-8.3'):
    source_files = [f[:-2] for f in os.listdir(object_dir) if f.endswith('.o')]
    
    # find all gcov data files
    for f in source_files:
        for ext in gcov_data_file_extensions:
            src = os.path.join(object_dir, f + ext)
            dst = os.path.join(object_dir, os.path.splitext(f)[0] + ext)
            if os.path.isfile(src):
                shutil.copy(src, dst)
    
    #call gcov on each source file
    cmd = [gcov_exe, '--object-directory', object_dir]
    relative_source_dir = os.path.relpath(source_dir)
    cmd.extend([os.path.join(relative_source_dir, f) for f in source_files])
    p = subprocess.Popen(cmd, universal_newlines=True,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.communicate()


def parse_arguments():
    p = argparse.ArgumentParser(
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('source_dir',
                   help="Directory with C/C++ source files.")
    p.add_argument('object_dir',
                   help="Directory with object files.")
    p.add_argument('--gcov', default='/pds/opt/gcc/bin/gcov-8.3',
                   help="gcov executable to use.")
    p.add_argument('--summary-file', default='summary.txt',
                   help=("Name of file that stores report summary "
                         "(will be placed in output directory, "
                         "see also other option)."))
    p.add_argument('-f', '--force', action='store_true',
                   help="Generate new report even if nothing has changed.")
    p.add_argument('-o', '--output-directory', default='coverage',
                   help="Directory to save gcov reports to.")
    p.add_argument('-v', action='store_true',
                   help="Show more detailed coverage information.")
    return p.parse_args()


if __name__ == '__main__':
    sys.exit(main())
