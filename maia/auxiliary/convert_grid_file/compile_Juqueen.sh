#!/bin/sh
mpicc -O2 convert_grid_file.c -o convert_grid_file -I/bgsys/local/parallel-netcdf/v1.5.0_bgclang/include /bgsys/local/parallel-netcdf/v1.5.0_bgclang/lib/libpnetcdf.a -lm 
