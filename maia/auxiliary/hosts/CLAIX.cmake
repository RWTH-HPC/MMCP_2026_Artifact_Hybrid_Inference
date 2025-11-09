################################################################################
# Set supported compilers
set(HOST_SUPPORTED_COMPILERS "GNU" "Intel")

################################################################################
# Set default compiler
set(HOST_DEFAULT_COMPILER "GNU")

################################################################################
# Set executable for each supported compiler (as 'HOST_COMPILER_<compiler>'))
set(HOST_COMPILER_EXE_Intel "mpiicpc")
set(HOST_COMPILER_EXE_GNU "mpicxx")

################################################################################
# Set default settings
# Include directories

set(INCLUDE_DIRS
)

# Library directories

set(LIBRARY_DIRS
)

################################################################################
# Set additional include/library directories or libraries for certain compilers

set(INCLUDE_DIRS_Intel
	    "~/libraries/fftw_intel/include"
)

set(LIBRARY_DIRS_Intel
	    "~/libraries/fftw_intel/lib"
)

# Library names
set(LIBRARY_NAMES
    "fftw3_mpi"
    "fftw3"
    "pnetcdf"
)

set(HOST_CONFIGURE_GNU
    "module load GCC/11.3.0 OpenMPI"
    "module load PnetCDF"
    "module load FFTW.MPI"
    "module load CMake"
)

set(HOST_CONFIGURE_Intel
    "module load intel/2022a impi"
    "module load PnetCDF"
    "module load CMake"
)

if(WITH_HDF5)
   set(HOST_CONFIGURE_INTEL
       "module load LIBRARIES"
       "module load hdf5/1.10.4"
       "module load intelmkl"
   )
endif()

# CANTERA
if(WITH_CANTERA)
  find_package(Threads REQUIRED)
  set(INCLUDE_DIRS ${INCLUDE_DIRS} ${SRC_DIR_ABS}/../include/canteraSoftLink/include)
  set(LIBRARY_DIRS ${LIBRARY_DIRS} ${SRC_DIR_ABS}/../include/canteraSoftLink/build/lib  )
  set(LIBRARY_NAMES ${LIBRARY_NAMES} "cantera")
  set(CANTERA_DATA ${CANTERA_DATA} ${SRC_DIR_ABS}/../include/canteraSoftLink/build/data)
endif()

################################################################################
# Set additional include/library directories or libraries for certain compilers

################################################################################
# Set (if necessary) configuration commands that should be executed when
# configuring MAIA.
set(HOST_PRESERVE_ENV "=")

