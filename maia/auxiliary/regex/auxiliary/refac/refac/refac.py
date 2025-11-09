#!/usr/bin/env python3

import argparse
import importlib
import os
import re
import shutil
import sys

def parse_arguments():
    # Get list of available modules (but exclude current script)
    scriptdir = os.path.dirname(os.path.realpath(__file__))
    modules = [m[:-3] for m in os.listdir(scriptdir) if m.endswith('.py') and
               m != os.path.basename(__file__)]

    # Set up parser
    p = argparse.ArgumentParser(
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('--in-place', '-i', action='store_true',
                   help="Modify files in-place.")
    p.add_argument('--stdout', '-1', action='store_true',
                   help="Print modified file to STDOUT "
                        "instead of <filename>_mod.")
    p.add_argument('module', help="Refactoring module to use.", choices=modules)
    p.add_argument('filename', nargs='*',
                   help="List of filenames to parse. If not specified, parse "
                        "all files in the current directory.")

    # Parse arguments
    args =  p.parse_args()

    # Sanity check for some arguments
    if args.in_place and args.stdout:
        sys.stderr.write("error: cannot use '--in-place' and '--stdout' at the "
                         "same time.\n")
        sys.exit(2)
    return args

# Computes the zfs file type for a given filename
def file_type(filename):
    filename = os.path.basename(filename)
    grid_files = ["zfscartesiangrid.h",
                  "zfscartesiangrid_hilbert.h",
                  "zfsgridgenpar.h",
                  "zfscartesiangrid.cpp",
                  "zfscartesiangrid_hilbert.cpp",
                  "zfscartesiangrid_inst_avg.cpp",
                  "zfscartesiangrid_inst_fv.cpp",
                  "zfscartesiangrid_inst_lb.cpp",
                  "zfscartesiangrid_inst_dg.cpp",
                  "zfsgridgenpar.cpp"]
    # Set type of cpp file
    is_header = False;
    is_cpp = False;
    is_source_file = True;
    if filename.endswith(".h"):
        is_header = True
    elif filename.endswith(".cpp"):
        is_cpp = True
    else:
        is_source_file = False
        return {}

    ftype = { 'grid':False, 'cell':False, 'block':False,
              'bnd':False, 'partcont':False, 'particle':False, 'lbm':False,
              'dg':False, 'fv':False, 'geo':False, 'interface':False,
              'avg':False, 'gridgenpar':False,
              'postproc':False, 'utility':False }

    for k,v in ftype.items():
        if k == 'grid':
            if filename in grid_files:
                ftype[k] = True
        elif k in filename:
            ftype[k] = True;

    if all(val == False for (key,val) in ftype.items()):
        ftype['utility'] = True
    return ftype

# extract delimited occurrence
def extract_first_delimited_occurrence(filestring, match):
    # if 'input_cells' in match['begin_delimiter']:
    #     print ("NEW CALL: input_cells matcher")

    # precondition: match must be a delimited match
    if not 'end_delimiter' in match or not 'begin_delimiter' in match:
        raise Exception("not a delimited match")

    # Find the delimiters:
    allowed_delimiters = [ ('[',']'), ('(',')') ]
    delimiters = None
    for b_delimiter,e_delimiter in allowed_delimiters:
        bd = b_delimiter if match['begin_delimiter'].endswith("\\"+b_delimiter) else None
        ed = e_delimiter if match['end_delimiter'].startswith("\\"+e_delimiter) else None
        if bd is None or ed is None:
            if not bd is ed:
                exit("mismatched delimiters: " + str(b_delimiter) + " " + str(e_delimiter))
        else: delimiters = (bd, ed)
    if delimiters is None: exit("delimiter not found")

    # Find first occurrence of the end delimiter
    pos_end_delimiter_b_m = re.search(match['end_delimiter'], filestring, re.MULTILINE)
    if pos_end_delimiter_b_m == None: return None

    pos_end_delimiter_b = pos_end_delimiter_b_m.start()
    pos_end_delimiter_e = pos_end_delimiter_b_m.end()
    pos_begin_delimiter_b = None
    pos_begin_delimiter_e = None


    iteration_no = 0
    while True:
        #if iteration_no > 10: exit("too many iterations")
        # if 'input_cells' in match['begin_delimiter']:
        #     print("it:" + str(iteration_no)
        #           + " | b_del_pos: [" + str(pos_begin_delimiter_b)
        #           + ", " + str(pos_begin_delimiter_e)  + "]"
        #           + " | e_del_pos: [" + str(pos_end_delimiter_b)
        #           + ", " + str(pos_end_delimiter_e)  + "]")
        #     if pos_begin_delimiter_b != None and pos_begin_delimiter_e != None:
        #         print("b_del_str: " + filestring[pos_begin_delimiter_b:pos_begin_delimiter_e])
        #     if pos_end_delimiter_b != None and pos_end_delimiter_e != None:
        #         print("e_del_str: " + filestring[pos_end_delimiter_b:pos_end_delimiter_e])

        # Find corresponding begin delimiter
        del_counter = 0
        # position of matching delimiter:
        pos_begin_delimiter_e = pos_end_delimiter_b
        for p in range(pos_end_delimiter_b, 0, -1):
            # if 'input_cells' in match['begin_delimiter']:
            #     print("dcount: " + str(del_counter) + " | p: " + str(p) + " | f[p]: " + filestring[p])
            if filestring[p] is delimiters[1]: # add 1 on end delimiter
                del_counter = del_counter + 1
            if filestring[p] is delimiters[0]: # substract 1 on begin delimiter
                del_counter = del_counter - 1
            if del_counter < 0: exit("wtf")
            if del_counter is 0:
                pos_begin_delimiter_e = p
                break
        if pos_begin_delimiter_e == None: exit("begin delimiter is none")
        if pos_begin_delimiter_e == pos_end_delimiter_b: exit("begin delimiter not found")

        # if 'input_cells' in match['begin_delimiter']:
        #     print("it:" + str(iteration_no)
        #           + " | b_del_pos: [" + str(pos_begin_delimiter_b)
        #           + ", " + str(pos_begin_delimiter_e) + "]"
        #           + " | e_del_pos: [" + str(pos_end_delimiter_b)
        #           + ", " + str(pos_end_delimiter_e)  + "]")
        #     if pos_begin_delimiter_b != None and pos_begin_delimiter_e != None:
        #         print("b_del_str: " + filestring[pos_begin_delimiter_b:pos_begin_delimiter_e])
        #     if pos_end_delimiter_b != None and pos_end_delimiter_e != None:
        #         print("e_del_str: " + filestring[pos_end_delimiter_b:pos_end_delimiter_e])

        # here pos_begin_delimiter_e points to the correct matching delimiter

        # Find if the begin matches the begin regular expression: Find all
        # occurrences of the regular expresion up to pos_begin_delimiter_e +1
        all_begin_occurrences = [(m.start(0), m.end(0)) for m in re.finditer(match['begin_delimiter'], filestring[0:pos_begin_delimiter_e+1], re.MULTILINE)]

        # if 'input_cells' in match['begin_delimiter']:
        #     print(str(all_begin_occurrences))

        if len(all_begin_occurrences) == 0:
            # if 'input_cells' in match['begin_delimiter']:
            #     print("NOT FOUND")
            # find next match with end delimiter
            offset = re.search(match['end_delimiter'], filestring[pos_end_delimiter_e:], re.MULTILINE)
            # not found (EOF) -> None
            if offset == None: return None

            pos_end_delimiter_b = pos_end_delimiter_e + offset.start()
            pos_end_delimiter_e = pos_end_delimiter_e + offset.end()
            # if 'input_cells' in match['begin_delimiter']:
            #     print("new pos_end_del: [" + str(pos_end_delimiter_b) + ", " + str(pos_end_delimiter_e))
            #     print("new end_del_str: " + filestring[pos_end_delimiter_b:pos_end_delimiter_e])
            iteration_no = iteration_no + 1

            continue

        pos_begin_delimiter_b_found = all_begin_occurrences[-1][0]
        pos_begin_delimiter_e_found = all_begin_occurrences[-1][1]

        # if 'input_cells' in match['begin_delimiter']:
        #     print(str(iteration_no)
        #           + " | b_del_f_pos: [" + str(pos_begin_delimiter_b_found)
        #           + ", " + str(pos_begin_delimiter_e_found) + "]")

        #     if pos_begin_delimiter_b_found != None and pos_begin_delimiter_e_found != None:
        #         print("found b_del_str: " + filestring[pos_begin_delimiter_b_found:pos_begin_delimiter_e_found])

        if pos_begin_delimiter_e_found == pos_begin_delimiter_e + 1:
            pos_begin_delimiter_b = pos_begin_delimiter_b_found
#             if 'input_cells' in match['begin_delimiter']:
#                 print("FOUND: it:" + str(iteration_no)
#                       + " | b_del_pos: [" + str(pos_begin_delimiter_b)
#                       + ", " + str(pos_begin_delimiter_e)  + "]"
#                       + " | e_del_pos: [" + str(pos_end_delimiter_b)
#                       + ", " + str(pos_end_delimiter_e)  + "]")
#                 if pos_begin_delimiter_b != None and pos_begin_delimiter_e != None:
#                     print("b_del_str: " + filestring[pos_begin_delimiter_b:pos_begin_delimiter_e])
#                 if pos_end_delimiter_b != None and pos_end_delimiter_e != None:
#                     print("e_del_str: " + filestring[pos_end_delimiter_b:pos_end_delimiter_e])
# co
            # found matches: return the match
            return {
                'occurrence_b':pos_begin_delimiter_b,
                'occurrence_e':pos_end_delimiter_e-1,
                'delimited_b':pos_begin_delimiter_e+1,
                'delimited_e':pos_end_delimiter_b,
                'occurrence_str':filestring[pos_begin_delimiter_b:pos_end_delimiter_e+1],
                'delimited_str':filestring[pos_begin_delimiter_e+1:pos_end_delimiter_b],
                'file_string':filestring,
            }
        else:
            # if 'input_cells' in match['begin_delimiter']:
            #     print("NOT FOUND")
            # find next match with end delimiter
            offset = re.search(match['end_delimiter'], filestring[pos_end_delimiter_e:], re.MULTILINE)
            # not found (EOF) -> None
            if offset == None: return None

            pos_end_delimiter_b = pos_end_delimiter_e + offset.start()
            pos_end_delimiter_e = pos_end_delimiter_e + offset.end()
            # if 'input_cells' in match['begin_delimiter']:
            #     print("new pos_end_del: [" + str(pos_end_delimiter_b) + ", " + str(pos_end_delimiter_e))
            #     print("new end_del_str: " + filestring[pos_end_delimiter_b:pos_end_delimiter_e])
            iteration_no = iteration_no + 1


def extract_first_token_occurrence(filestring, match):
    # Check if a correct matching object was passed
    if not 'token' in match:
        raise Exception("not a token match")

    # Find first occurrence of the token
    m = re.search(match['token'], filestring, re.MULTILINE)

    # Return None if token was not found
    if m is None:
        return None

    # Return proper occurrence dict
    return {
        'occurrence_b': m.start(),
        'occurrence_e': m.end()-1,
        'delimited_b': None,
        'delimited_e': None,
        'occurrence_str': m.group(),
        'delimited_str': None,
        'file_string': filestring,
    }


# extracts the first occurrence inside filestring that matches match
def first_occurrence(filestring, match):
    # right now only the delimited version is implemented
    if 'end_delimiter' in match and 'begin_delimiter' in match:
        oc = extract_first_delimited_occurrence(filestring, match)
    elif 'token' in match:
        oc = extract_first_token_occurrence(filestring, match)
    else:
        exit("not implemented")

    # if oc is None:
    #     print("no occurrence found for matcher: " + str(match))
    #else:
        # print("occurrence found:")
        # for k,v in oc.items():
        #     print(k + " : " + str(v))
    return oc

# rewrite the first occurrence in a filestring that matches match
def rewritte_occurrence(occurrence, match, ftype):
    before_occurrence = occurrence['file_string'][:occurrence['occurrence_b']]
    after_occurrence = occurrence['file_string'][occurrence['occurrence_e']+1:]
    rewritten_occurrence = match['rewrite_function'](occurrence, match, ftype)
    return before_occurrence + rewritten_occurrence + after_occurrence


# For every file on the directory, loads the file into a single string. It then
# iterates over the occurrences of the "match" pattern within the string
# replacing them according to the "rewritte_occurrence" function.
def main():
    # Get command line arguments
    args = parse_arguments()

    # Get correct module
    module = importlib.import_module(args.module)

    # Get list of filenames
    if len(args.filename) == 0:
        filenames = [f for f in os.listdir(".")]
    else:
        filenames = args.filename
        not_files = [not os.path.isfile(f) for f in filenames]
        if any(not_files):
            err = "error: the following arguments are no files:\n"
            for i,s in enumerate(not_files):
                if not s:
                    continue
                err += "  - " + filenames[i] + '\n'
            sys.stderr.write(err)
            sys.exit(2)

    # For each file in the current directory:
    for filename in filenames:
        # Compute the file type:
        ftype = file_type(filename)
        # Skip non-source files
        if len(ftype) == 0:
            print("[SKIP] " + filename + " (not a ZFS source code file)")
            continue
        else:
            print("["
                  + "".join([ k + " " for k,v in ftype.items() if v]).rstrip()
                  + "] " + filename)

        # Load the whole file into a single string
        file_string = open(filename, 'r', newline='').read()

        # Set the output stream (either stdout or a file stream)
        tmpfile = filename + "_mod"
        output_stream = sys.stdout if args.stdout else open(tmpfile, 'w')

        # Remember if a file was changed at all to avoid unnecessary copying
        file_has_changed = False

        # Process each occurrence into the new file string by appending the
        # already processed occurrences with the current one to allow
        # backtracking over previous occurrences (it is performed for every
        # matcher independently):
        for match in module.matchers:
            occurrence = first_occurrence(file_string, match)
            while occurrence != None:
                file_has_changed = True
                file_string = rewritte_occurrence(occurrence, match, ftype)
                occurrence = first_occurrence(file_string, match)

        # Process each token
        for token in module.tokens:
            # Skip if wrong file type
            has_correct_filetype = False
            for ft, status in ftype.items():
                if status and ft in token['ftype']:
                    has_correct_filetype = True
            if not has_correct_filetype:
                continue

            # Replace tokens
            file_string, no_subs = re.subn(token['pattern'],
                                           token['replacement'],
                                           file_string,
                                           re.MULTILINE)
            if no_subs > 0:
                file_has_changed = True

        # Writes the new file string to the selected output stream:
        output_stream.write(file_string)

        # Replace file with _mod file if it has changed
        if args.in_place:
            # Move modified file to original file if there were changes,
            # otherwise just delete the modification file
            if file_has_changed:
                shutil.move(tmpfile, filename)
            else:
                os.remove(tmpfile)

if __name__ == '__main__':
    sys.exit(main())
