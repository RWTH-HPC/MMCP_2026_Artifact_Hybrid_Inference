################################################################################
# Set supported compilers
set(HOST_SUPPORTED_COMPILERS "GNU" "AMD" "Intel")

################################################################################
# Set default compiler
set(HOST_DEFAULT_COMPILER "GNU")

################################################################################
# Set executable for each supported compiler (as 'HOST_COMPILER_<compiler>'))
set(HOST_COMPILER_EXE_GNU   "mpicxx")
set(HOST_COMPILER_EXE_AMD   "mpicxx")
set(HOST_COMPILER_EXE_Intel "mpicxx")
set(HOST_COMPILER_EXE_PGI   "mpicxx")

################################################################################
# Set default settings
# Include directories
set(INCLUDE_DIRS)

# Library directories
set(LIBRARY_DIRS)

# Library names
set(LIBRARY_NAMES
  "fftw3_mpi"
  "fftw3"
  "pnetcdf"
)

# CANTERA
if(WITH_CANTERA)
  find_package(Threads REQUIRED)
  set(INCLUDE_DIRS ${INCLUDE_DIRS} ${SRC_DIR_ABS}/../include/canteraSoftLink/include)
  set(LIBRARY_DIRS ${LIBRARY_DIRS} ${SRC_DIR_ABS}/../include/canteraSoftLink/build/lib  )
  set(LIBRARY_NAMES ${LIBRARY_NAMES} "cantera")
  set(CANTERA_DATA ${CANTERA_DATA} ${SRC_DIR_ABS}/../include/canteraSoftLink/build/data)
endif()

set(LIBRARY_DIRS_Intel
  # TODO UPDATE
  "/opt/hlrs/spack/rev-003_2020-03-03/parallel-netcdf/1.12.0-intel-19.1.0-bax47cmw/lib"
  "/opt/hlrs/spack/rev-003_2020-03-03/fftw/3.3.8-intel-19.1.0-p2vkuugf/lib"
)

set(LIBRARY_DIRS_AMD
  # TODO UPDATE
  "/opt/hlrs/spack/rev-003_2020-03-03/parallel-netcdf/1.12.0-clang-9.0.0-asejiopq/lib"
  "/opt/hlrs/spack/rev-003_2020-03-03/fftw/3.3.8-clang-9.0.0-3ysve2s3/lib"
)

################################################################################
# Set (if necessary) configuration commands that should be executed when
# configuring MAIA
set(HOST_CONFIGURE
  "module unload gcc"
  "module load python"
)

set(HOST_CONFIGURE_GNU
  "module load gcc/10.2.0"
  "module load parallel-netcdf/1.12.1"
  "module load fftw/3.3.8"
  "module load hdf5"
)

# TODO Untested!
set(HOST_CONFIGURE_AMD
  "module load aocc/2.1.0"
  "module load parallel-netcdf/1.12.0"
  "module load hdf5"
)

# TODO Untested!
set(HOST_CONFIGURE_Intel
  "module load intel/19.1.3"
  "module load parallel-netcdf/1.12.0"
  "module load fftw/3.3.8"
  "module load hdf5"
)

if(WITH_HDF5)
  set(LIBRARY_NAMES ${LIBRARY_NAMES} "hdf5")
  set(LIBRARY_NAMES ${LIBRARY_NAMES} "z")
  set(LIBRARY_NAMES ${LIBRARY_NAMES} "dl")
endif()

################################################################################
# Specify environment variables that should be preserved from configure time
# If you add '=', all environment variables will be preserved
set(HOST_DEFAULT_OPTIONS "WITH_HDF5")
set(HOST_PRESERVE_ENV "=")
