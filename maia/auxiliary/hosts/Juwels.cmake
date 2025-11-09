################################################################################
# Set supported compilers
set(HOST_SUPPORTED_COMPILERS "GNU" "Intel" "NVHPC")

################################################################################
# Set default compiler
set(HOST_DEFAULT_COMPILER "GNU")

################################################################################
# Set executable for each supported compiler (as 'HOST_COMPILER_<compiler>'))
set(HOST_COMPILER_EXE_GNU "mpicxx")
set(HOST_COMPILER_EXE_Intel "mpicxx")
set(HOST_COMPILER_EXE_NVHPC "mpicxx")

################################################################################
# Set default settings
# Include directories
set(INCLUDE_DIRS)

# Library directories
set(LIBRARY_DIRS)

# Library names
set(LIBRARY_NAMES
    "cudart"
    "fftw3_mpi"
    "fftw3"
    "pnetcdf"
    "mpi"
    )

#PSTL - cuda version
set(CUDA_VERSION "11.4")

################################################################################
# Set additional include/library directories or libraries for certain compilers

set(LIBRARY_DIRS_GNU
  "$ENV{EBROOTFFTW}/lib"
  "$ENV{EBROOTPARALLELMINNETCDF}/lib"  
  "$ENV{EBROOTPSCOM}/lib"
  "$ENV{EBROOTHDF5}/lib"
)

set(LIBRARY_DIRS_Intel
  "$ENV{EBROOTFFTW}/lib"
  "$ENV{EBROOTPARALLELMINNETCDF}/lib"
  "$ENV{EBROOTPSCOM}/lib"
  "$ENV{EBROOTHDF5}/lib"
)

set(LIBRARY_DIRS_NVHPC
  "$ENV{EBROOTPARALLELMINNETCDF}/lib"
  "$ENV{EBROOTOPENMPI}/lib"
  "$ENV{EBROOTNVHPC}/Linux_x86_64/21.9/cuda/11.4/lib64"
)
set(INCLUDE_DIRS_NVHPC
  "$ENV{EBROOTNVHPC}/Linux_x86_64/21.9/cuda/11.4/include"
)
set(HOST_CONFIGURE_NVHPC
  "module load Stages/2020 NVHPC/21.9-GCC-10.3.0 OpenMPI/4.1.1 parallel-netcdf Boost.Python FFTW/3.3.8"
)
#PSTL - cuda version
set(CUDA_VERSION "11.4")


################################################################################
# Set (if necessary) configuration commands that should be executed when
# configuring MAIA.
set(HOST_CONFIGURE_GNU
  "module load GCC ParaStationMPI FFTW parallel-netcdf Python HDF5"
)
set(HOST_CONFIGURE_Intel
  "module load Intel ParaStationMPI FFTW parallel-netcdf intel-para HDF5"
)
set(HOST_CONFIGURE_NVHPC
  "module load Stages/2020 NVHPC/21.9-GCC-10.3.0 OpenMPI/4.1.1 FFTW/3.3.8 parallel-netcdf Boost.Python"
)
################################################################################
# Specify environment variables that should be preserved from configure time
# If you add '=', all environment variables will be preserved

