################################################################################
# Set supported compilers
# set(HOST_SUPPORTED_COMPILERS "GNU" "Clang" "Intel" "PGI" "GNU_nightly")
set(HOST_SUPPORTED_COMPILERS "PGI")

################################################################################
# Set default compiler
set(HOST_DEFAULT_COMPILER "PGI")

################################################################################
# Set executable for each supported compiler (as 'HOST_COMPILER_<compiler>'))
# set(HOST_COMPILER_EXE_GNU
#     "/pds/opt/gcc-4.7-2/bin/g++")
# set(HOST_COMPILER_EXE_GNU_nightly
#     "/home/mic/.pool/gcc_nightly/bin/g++")
# set(HOST_COMPILER_EXE_Clang
#     "/home/mic/.pool/llvm_nightly/bin/clang++")
# set(HOST_COMPILER_EXE_Intel
#     "/home/mic/.pool/intel/14.0.2/compiler/bin/icpc")
set(HOST_COMPILER_EXE_PGI
     "$ENV{PGI_HOME}/bin/pgc++")

################################################################################
# Set default settings
# Include directories
set(INCLUDE_DIRS
    "$ENV{FFTLIB}/include"
    "$ENV{PARALLEL_NETCDF_HOME}/include"
    "$ENV{MPI_HOME}/include"
)

# Library directories
set(LIBRARY_DIRS
    "$ENV{FFTLIB}/lib"
    "$ENV{PARALLEL_NETCDF_HOME}/lib"
    "$ENV{MPI_HOME}/lib"
    "$ENV{LIBNUMA_HOME}/lib64"
)

# Library names
set(LIBRARY_NAMES
    "fftw3"
    "pnetcdf"
    "mpi_cxx"
    "mpi"
    "rt"
    "curl"
    "numa"
)

################################################################################
# Set additional include/library directories or libraries for certain compilers
# Additional info for clang's libcxx
# set(INCLUDE_DIRS_Clang
#     "/home/mic/.pool/libcxx_nightly/include/c++/v1"
# )
# set(LIBRARY_DIRS_Clang
#     "/home/mic/.pool/libcxx_nightly/lib"
# )
# set(LIBRARY_NAMES_Clang
#     "c++"
# )

#ZLIB
# if(WITH_ZLIB)
#   set(INCLUDE_DIRS ${INCLUDE_DIRS} "/home/mic/.pool/zlib-1.2.8/include")
#   set(LIBRARY_DIRS ${LIBRARY_DIRS} "/home/mic/.pool/zlib-1.2.8/lib")
#   set(LIBRARY_NAMES ${LIBRARY_NAMES} "z")
# endif()

################################################################################
# Enable clang static analyzer for this host
# set(CLANG_STATIC_ANALYZER_DIR
#     "/home/mic/.pool/.src/llvm/tools/clang/tools/scan-build")

################################################################################
# Set (if necessary) configuration commands that should be executed when
# configuring MAIA.
# set(HOST_CONFIGURE_GNU
#     "export LD_LIBRARY_PATH=/home/mic/.pool/gcc-4.8.2/lib64:$LD_LIBRARY_PATH")
# set(HOST_CONFIGURE_GNU_nightly
#     "export LD_LIBRARY_PATH=/home/mic/.pool/gcc_nightly/lib64:$LD_LIBRARY_PATH")
# set(HOST_CONFIGURE_Intel
#     "export INTEL_LICENSE_FILE=50017@license2.rz.rwth-aachen.de")
# set(HOST_CONFIGURE_PGI
#     "export PGROUPD_LICENSE_FILE=50000@license2.rz.rwth-aachen.de")

################################################################################
# Specify environment variables that should be preserved from configure time
# If you add '=', all environment variables will be preserved
# set(HOST_PRESERVE_ENV_Intel
#     "INTEL_LICENSE_FILE"
# )
# set(HOST_PRESERVE_ENV_PGI
#     "PGROUPD_LICENSE_FILE"
# )
