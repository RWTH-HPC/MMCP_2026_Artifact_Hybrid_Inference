#!/bin/sh
module load PrgEnv-gnu
module load cray-parallel-netcdf
cc -std=c99 -O2 convert_grid_file.c -o convert_grid_file -lm
