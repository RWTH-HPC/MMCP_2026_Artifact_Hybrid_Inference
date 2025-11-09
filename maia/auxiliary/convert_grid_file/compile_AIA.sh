#!/bin/sh
mpicc -std=c99 convert_grid_file.c -o convert_grid_file -I/pds/opt/parallel-netcdf/include /pds/opt/parallel-netcdf/lib/libpnetcdf.a -lm
