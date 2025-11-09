################################################################################
# Set supported compilers
set(HOST_SUPPORTED_COMPILERS "GNU" "Intel")

################################################################################
# Set default compiler
set(HOST_DEFAULT_COMPILER "GNU")

################################################################################
# Set executable for each supported compiler (as 'HOST_COMPILER_<compiler>'))
set(HOST_COMPILER_EXE_GNU "mpicxx")
set(HOST_COMPILER_EXE_Intel "mpicxx")

################################################################################
# Set default settings

# Library names
if(WITH_KNLBOOSTER)
 set(LIBRARY_NAMES
    "fftw3_mpi"
    "fftw3"
    "pnetcdf"
    )
else()
set(LIBRARY_NAMES
    "fftw3_mpi"
    "fftw3"
    "pnetcdf"
    "pscom"
    )
endif() 

if(WITH_HDF5)
   set(LIBRARY_NAMES ${LIBRARY_NAMES} "dl")
   set(LIBRARY_NAMES ${LIBRARY_NAMES} "hdf5")
   set(LIBRARY_NAMES ${LIBRARY_NAMES} "hdf5_hl")
   set(LIBRARY_NAMES ${LIBRARY_NAMES} "sz")
   set(LIBRARY_NAMES ${LIBRARY_NAMES} "z")
endif()

################################################################################
# Set additional include/library directories or libraries for certain compilers

set(LIBRARY_DIRS_GNU
  "$ENV{EBROOTFFTW}/lib"
  "$ENV{EBROOTPARALLELMINNETCDF}/lib"
  "$ENV{EBROOTPSCOM}/lib"
)

if(WITH_HDF5IOLIB)
   set(LIBRARY_DIRS_GNU
     ${LIBRARY_DIRS_GNU}
     "$ENV{EBROOTHDF5}/lib"
     "$ENV{EBROOTSZIP}/lib"
     "$ENV{EBROOTZLIB}/lib"
   )
endif()


set(LIBRARY_DIRS_Intel
  "$ENV{EBROOTFFTW}/lib"
  "$ENV{EBROOTPARALLELMINNETCDF}/lib"
  "$ENV{EBROOTPSCOM}/lib"
)

################################################################################
# Set (if necessary) configuration commands that should be executed when
# configuring MAIA.
set(HOST_CONFIGURE_GNU
  "module load GCC ParaStationMPI FFTW parallel-netcdf HDF5"
)
set(HOST_CONFIGURE_Intel
  "module load Intel ParaStationMPI FFTW parallel-netcdf intel-para"
)

################################################################################
# Specify environment variables that should be preserved from configure time
# If you add '=', all environment variables will be preserved
set(HOST_PRESERVE_ENV "=")
