################################################################################
# Set supported compilers
set(HOST_SUPPORTED_COMPILERS "GNU" "Cray" "Intel" "PGI")

################################################################################
# Set default compiler
set(HOST_DEFAULT_COMPILER "GNU")

################################################################################
# Set executable for each supported compiler (as 'HOST_COMPILER_<compiler>'))
set(HOST_COMPILER_EXE_GNU   "CC")
set(HOST_COMPILER_EXE_Cray  "CC")
set(HOST_COMPILER_EXE_Intel "CC")
set(HOST_COMPILER_EXE_PGI   "CC -host=x86_64-suse-linux-gnu")

################################################################################
# Set default settings
# Include directories
set(INCLUDE_DIRS)

# Library directories
set(LIBRARY_DIRS)

# Library names
set(LIBRARY_NAMES)

################################################################################
# Set (if necessary) configuration commands that should be executed when
# configuring MAIA
set(HOST_CONFIGURE
    "module unload PrgEnv-intel PrgEnv-cray PrgEnv-pgi PrgEnv-gnu"
    "module unload fftw cray-netcdf cray-netcdf-hdf5parallel"
    "module unload cray-hdf5-parallel cray-parallel-netcdf"
)

set(HOST_CONFIGURE_GNU
    "module load PrgEnv-gnu"
    "module load cray-fftw cray-parallel-netcdf cray-hdf5-parallel"
)

set(HOST_CONFIGURE_Intel
    "module load PrgEnv-intel"
    "module load cray-fftw cray-parallel-netcdf cray-hdf5-parallel"
)

set(HOST_CONFIGURE_Cray
    "module load PrgEnv-cray"
    "module load gcc"
    "module load cray-fftw cray-parallel-netcdf cray-hdf5-parallel"
)


#HDF5IOLIB
if(WITH_HDF5IOLIB)
  set(INCLUDE_DIRS ${INCLUDE_DIRS} "/users/malbers/local_opt/iolib/build")
  set(LIBRARY_DIRS ${LIBRARY_DIRS} "/users/malbers/local_opt/iolib/build")
  set(LIBRARY_NAMES ${LIBRARY_NAMES} "maiaio")
  set(LIBRARY_NAMES ${LIBRARY_NAMES} "dl")
endif()

################################################################################
# Specify environment variables that should be preserved from configure time
# If you add '=', all environment variables will be preserved
set(HOST_PRESERVE_ENV "=")
