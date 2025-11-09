zfs_refactor
============

refactoring tool for zfs

Usage
-----

Go to the 'src/' directory of ZFS, then execute

    path/to/zfs_refactor/zfsrefactoring.py module [file [file]...]

The `module` must be the name of the regex module you would like to use. If no
`file` is specified, all files in the current directory are used. By default,
all replacements go into a file called `<original_name>_mod`. If you want to
print to the shell instead, use the `--stdout`/`-1` flag. If you want to
change the file in-place, use the `--in-place`/`-i` flag. To see a help
message, just call it with the `--help`/`-h` flag.
