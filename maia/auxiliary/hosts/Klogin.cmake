################################################################################
# Set supported compilers
set(HOST_SUPPORTED_COMPILERS "Fujitsu")

################################################################################
# Set default compiler
set(HOST_DEFAULT_COMPILER "Fujitsu")

################################################################################
# Set executable for each supported compiler (as 'HOST_COMPILER_<compiler>'))
set(HOST_COMPILER_EXE_Fujitsu "mpiFCCpx")

################################################################################
# Set default settings
# Include directories
set(INCLUDE_DIRS_Fujitsu
    "$ENV{FFTW_HOME}/include"
    "$ENV{PNETCDF_HOME}/include"
)

# Library directories
set(LIBRARY_DIRS_Fujitsu
    "$ENV{FFTW_HOME}/lib64"
    "$ENV{PNETCDF_HOME}/lib"
)

# Library names
set(LIBRARY_NAMES_Fujitsu
    "fftw3_mpi"
    "fftw3"
    "pnetcdf"
)

set(HOST_CONFIGURE_Fujitsu
    "export LC_ALL=en_US"
    "export LANG=C"
)

set(HOST_CONFIGURE_Fujitsu
    "source /home/system/Env_base"
)

################################################################################
# Specify environment variables that should be preserved from configure time
# If you add '=', all environment variables will be preserved
set(HOST_PRESERVE_ENV "=")
