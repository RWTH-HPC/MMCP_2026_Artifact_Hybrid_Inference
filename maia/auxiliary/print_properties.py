#!/usr/bin/env python3
import argparse, collections, os, re, subprocess, sys
documented_properties = None
zfs_properties = None

def cpp_src_files(path):
    return [os.path.join(dpath, f) for (dpath, _, fnames) in os.walk(path) \
            for f in fnames if os.path.splitext(f)[1] in [".cpp", ".c", ".h", ".hpp"]] 

def properties_in_file(filename):
    with open(filename, 'r') as f: data = f.read().replace('\n', ' ').replace(' ', '')
    return collections.Counter(re.findall('ZFSContext::getProperty\(\s*"([0-9a-zA-Z\_\-]+)"', data))

def is_zfs_root_dir(path):
    for d in ['auxiliary', 'doc', 'src', '.svn']:
        if not d in next(os.walk(path))[1]:
            print "The path needs to point to the ZFS root directory (e.g. `trunk/`).\n The current path is: %s" %path
            sys.exit(1)

def is_testcases_root_dir(path):
    for d in ['FV', 'DG', 'LBM', 'GRID']:
        if not d in next(os.walk(path))[1]:
            print "The testcase path needs to point to the testcases root directory (e.g. `testcases/`).\n The current path is: %s" %path
            sys.exit(1)

def find_cdl_files(path):
    fs = set()
    for root, subdir, files in os.walk(path):
        for f in files:
            filename, file_extension = os.path.splitext(f)
            if '.cdl' in file_extension and not 'geometry' in filename:
                a = os.path.abspath(os.path.join(root, f))
                if a not in files:
                    fs.add(a)
    return fs

def find_used_properties_(f):
    with open(f, 'r') as f: data = f.readlines()
    variable_lines = []
    reading_variables = False
    
    for l in data:
        if reading_variables is False and "variables:" in l:
            reading_variables = True
            continue
        if reading_variables is True and "data:" in l:
            break
        if reading_variables is False:
            continue
        variable_lines.append(l.split('//')[0])
    variable_lines = ' '.join(variable_lines)
    variable_lines = variable_lines.replace('int', '').replace('double', '').replace('char', '').replace(',', ' ').replace(';', ' ').split(' ')
    variable_lines = list(set([l.strip().split('(')[0].split('.')[0] for l in variable_lines if len(l.strip()) > 0]))
    return variable_lines
    
def find_used_properties(cdl_files):
    m = {}
    for f in cdl_files:
        used_properties = find_used_properties_(f)
        m[f] = used_properties
    return m

def find_unused_properties(zfs_properties, cdl_files):
    used_properties = find_used_properties(cdl_files)
    unused_properties = {}
    for (property_file,properties) in used_properties.iteritems():
        for p in properties:
            if p in zfs_properties:
                continue
            if property_file not in unused_properties:
                l = []
                l.append(p)
                unused_properties[property_file] = l
            else:
                unused_properties[property_file].append(p)
    return unused_properties
                                
def stats(properties):
    total = len(properties)
    if documented_properties:
        doc = len([x for x,c in properties.iteritems() if documented_properties[x] == c and zfs_properties[x] == c])
        undoc = len([x for x,_ in properties.iteritems() if x not in documented_properties])
        maybe = total - doc - undoc
        return "total: %s, documented: %s, undocumented: %s, maybe undocumented: %s" \
            %(str(total), str(doc), str(undoc), str(maybe))
    else:
        return "total: %s" % str(total)

def print_property(x, c):
    st = ""            
    if documented_properties:
        if documented_properties[x] == c:
            st = ""
        elif not x in documented_properties:
            st = "undocumented"
        else:
            st = "maybe undocumented [current: %s, zfs: %s, docs: %s] " \
                 %(str(c), str(zfs_properties[x]), str(documented_properties[x]))
    if len(st) > 0: st = "  (%s)" %st 
    print "  %s%s" %(x, st)

def print_properties(properties):
        for k,v in properties.iteritems(): print_property(k,v)
        
def print_undoc_properties(properties):
        numberUndoc = 0
        for k,v in properties.iteritems():
          if k not in documented_properties and k + 'N' not in documented_properties:
            print_property(k,v)
            numberUndoc += 1
        print "Number:", str(numberUndoc)
  
def build_zfs_documentation(src_dir):
    doc = subprocess.Popen(['cd %s && ./configure.py 1 1 && make doc' % src_dir],
                           stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    lines = doc.stdout.readlines() # TODO: doc.returncode is always !== 0 due to doxygen warnings...
    
def load_documented_properties(src_dir, build_docs):
    print "Scanning documentation for documented properties..."
    if build_docs: build_zfs_documentation(src_dir)
    html_dir = os.path.join(src_dir, 'doc/code/doxygen/html')
    if not os.path.exists(html_dir):
        print "The Doxygen documentation is not built. Did you forget to `make doc` or use --build-docs?"
        sys.exit(1)
    else:
        print "Doxygen documentation found."
    property_file_name = None
    for fname in next(os.walk(html_dir))[2]:
        # property_page1 / PropertyPage1 / etc.
        if 'roperty' in fname and 'age' in fname:
            property_file_name = fname
            break
    if property_file_name is None:
        print "Cannot find the Doxygen property file. Maybe you are using an unsupported version of Doxygen?"
        sys.exit(1)
    print property_file_name
    property_file = os.path.join(html_dir, property_file_name)
    with open(property_file, 'r') as f: data = f.read().replace('\n', ' ')
    global documented_properties
    documented_properties = collections.Counter(re.findall('<a class="anchor" id="([0-9a-zA-Z\_\-]+)">', data))

def main():
    parser = argparse.ArgumentParser(description="Prints ZFS properties")
    parser.add_argument("zfs_root_path", help="path to zfs source code directory (e.g. trunk/)", type=str)
    parser.add_argument("--all", help="prints all properties", action='store_true', dest="all")
    parser.add_argument("--undoc", help="prints all undocumented properties", action='store_true', dest="undoc")
    parser.add_argument("--per-file", help="prints all properties sorted per file", action='store_true', dest="per_file")
    parser.add_argument("--scan-docs", help="scans ZFS documentation (annotates undocumented properties)", action='store_true', dest="scan_docs")
    parser.add_argument("--build-docs", help="builds ZFS documentation", action='store_true', dest="build_docs", default=False)
    parser.add_argument("--stats",   help="prints property statistics", action='store_true', dest="stats")
    parser.add_argument("--omit-empty-files", help="omits files without properties", action='store_true', dest="omit_empty")
    parser.add_argument("--unused", help="prints properties that are not used in the testcases", dest="testcases", type=str,
                        metavar='PATH_TO_TESTCASES')
    args = parser.parse_args()

    is_zfs_root_dir(args.zfs_root_path)

    if not args.all and not args.per_file and not args.stats and not args.testcases and not args.undoc:
        print "ERROR: Nothing todo. Use either --all, --per-file, --stats, --undoc or --unused."
        sys.exit(0)

    file_properties =  [(f, properties_in_file(f)) for f in cpp_src_files(os.path.join(args.zfs_root_path, 'src/'))]

    if args.scan_docs or args.undoc: load_documented_properties(args.zfs_root_path, args.build_docs)
    global zfs_properties
    zfs_properties = sum([vs for _,vs in file_properties], collections.Counter()) 
    
    if args.undoc:
        print "UNDOCUMENTED PROPERTIES:"
        print_undoc_properties(zfs_properties)
    
    
    if args.all:
        print "ALL PROPERTIES (%s):" % stats(zfs_properties)
        print_properties(zfs_properties)

            
    if args.per_file:
        for k,v in file_properties:
            if args.omit_empty and len(v) == 0: continue
            print "%s (%s):" % (str(k), stats(v))
            print_properties(v)
                        
    if args.stats: print "STATS: %s" % stats(zfs_properties)

    if args.testcases:
        is_testcases_root_dir(args.testcases)
        cdl_files = find_cdl_files(args.testcases)
        unused_properties = find_unused_properties(zfs_properties, cdl_files)
        for (k,v) in unused_properties.iteritems():
            print "%s:" %k
            for i in v:
                print "    %s" %i

if __name__ == '__main__':  main()
