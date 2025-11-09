################################################################################
# Set supported compilers
set(HOST_SUPPORTED_COMPILERS "GNU")

################################################################################
# Set default compiler
set(HOST_DEFAULT_COMPILER "GNU")

################################################################################
# Set executable for each supported compiler (as 'HOST_COMPILER_<compiler>'))
set(HOST_COMPILER_EXE_GNU       "mpicxx")

################################################################################
# Set default settings
# Include directories
set(INCLUDE_DIRS
    "$ENV{HOME}/libraries/parallel-netcdf/include"
    "$ENV{HOME}/libraries/fftw/include"
)

# Library directories
set(LIBRARY_DIRS
    "$ENV{HOME}/libraries/parallel-netcdf/lib"
    "$ENV{HOME}/libraries/fftw/lib"
)

# Library names
set(LIBRARY_NAMES
    "fftw3"
    "pnetcdf"
    "pthread"
    "rt"
)

################################################################################
# Set additional include/library directories or libraries for certain compilers

################################################################################
# Set (if necessary) configuration commands that should be executed when
# configuring MAIA.
