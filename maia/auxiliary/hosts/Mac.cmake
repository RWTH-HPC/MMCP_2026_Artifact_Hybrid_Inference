################################################################################
# Set supported compilers
set(HOST_SUPPORTED_COMPILERS "GNU" "Clang")

################################################################################
# Set default compiler
set(HOST_DEFAULT_COMPILER "GNU")

################################################################################
# Set executable for each supported compiler (as 'HOST_COMPILER_<compiler>'))
set(HOST_COMPILER_EXE_GNU
    "/pds/opt/gcc/bin/g++")
set(HOST_COMPILER_EXE_Clang
    "mpic++")

################################################################################
# Set default settings
# Include directories
set(INCLUDE_DIRS
    "/pds/opt/netcdf/include"
    "/pds/opt/fftw/include"
    "/pds/opt/parallel-netcdf/include"
    "/pds/opt/openmpi/include"
)

# Library directories
set(LIBRARY_DIRS
    "/pds/opt/netcdf/lib"
    "/pds/opt/fftw/lib"
    "/pds/opt/parallel-netcdf/lib"
    "/pds/opt/openmpi/lib"
)

# Library names
set(LIBRARY_NAMES
    "netcdf"
    "fftw3"
    "pnetcdf"
    "mpi_cxx"
    "mpi"
    "curl"
)

################################################################################
# Set additional include/library directories or libraries for certain compilers
# Additional info for clang's libcxx
#set(INCLUDE_DIRS_Clang
#    "/home/mic/.pool/libcxx_nightly/include/c++/v1"
#)
#set(LIBRARY_DIRS_Clang
#    "/home/mic/.pool/libcxx_nightly/lib"
#)
#set(LIBRARY_NAMES_Clang
#    "c++"
#)

#ZLIB
#if(WITH_ZLIB)
#  set(INCLUDE_DIRS ${INCLUDE_DIRS} "/home/mic/.pool/zlib-1.2.8/include")
#  set(LIBRARY_DIRS ${LIBRARY_DIRS} "/home/mic/.pool/zlib-1.2.8/lib")
#  set(LIBRARY_NAMES ${LIBRARY_NAMES} "z")
#endif()

#HDF5IOLIB
if(WITH_HDF5IOLIB)
  set(INCLUDE_DIRS ${INCLUDE_DIRS} "/pds/opt/hdf5/lib")
  set(INCLUDE_DIRS ${INCLUDE_DIRS} "/home/marian/Documents/Masterarbeit/iolib")
  set(LIBRARY_DIRS ${LIBRARY_DIRS} "/pds/opt/hdf5/lib")
  set(LIBRARY_DIRS ${LIBRARY_DIRS} "/home/marian/Documents/Masterarbeit/iolib")
  set(LIBRARY_NAMES ${LIBRARY_NAMES} "maiaio")
  set(LIBRARY_NAMES ${LIBRARY_NAMES} "hdf5_hl")
  set(LIBRARY_NAMES ${LIBRARY_NAMES} "hdf5")
  set(LIBRARY_NAMES ${LIBRARY_NAMES} "z")
  set(LIBRARY_NAMES ${LIBRARY_NAMES} "dl")	
endif()

################################################################################
# Enable clang static analyzer for this host
#set(CLANG_STATIC_ANALYZER_DIR
#    "/home/mic/.pool/.src/llvm/tools/clang/tools/scan-build")

################################################################################
# Set (if necessary) configuration commands that should be executed when
# configuring MAIA.
#set(HOST_CONFIGURE_GNU
#    "export LD_LIBRARY_PATH=/home/mic/.pool/gcc/4.9.1/lib64:$LD_LIBRARY_PATH")

