################################################################################
# Set supported compilers
set(HOST_SUPPORTED_COMPILERS "GNU" "Clang" "Intel" "PGI" "NVHPC")

################################################################################
# Set default compiler
set(HOST_DEFAULT_COMPILER "GNU")

################################################################################
# Set executable for each supported compiler (as 'HOST_COMPILER_<compiler>'))
set(HOST_COMPILER_EXE_GNU
    "$ENV{MAIA_COMPILER_HOME}/gcc-11.1.0/bin/g++")
set(HOST_COMPILER_EXE_Clang
    "$ENV{MAIA_COMPILER_HOME}/llvm-9.0.1/bin/clang++")
set(HOST_COMPILER_EXE_Intel
    "$ENV{MAIA_COMPILER_HOME}/intel/oneapi/compiler/latest/linux/bin/icpx")
set(HOST_COMPILER_EXE_PGI
    "$ENV{MAIA_COMPILER_HOME}/pgi-19.10/bin/pgc++")
set(HOST_COMPILER_EXE_NVHPC
    "$ENV{MAIA_COMPILER_HOME}/nvhpc-23.7/Linux_x86_64/23.7/compilers/bin/nvc++")

set(HOST_DOXYGEN_EXE
    "$ENV{MAIA_COMPILER_HOME}/bin/doxygen")

################################################################################
# Set default settings
# Include directories
set(INCLUDE_DIRS
    "$ENV{MAIA_LIB_HOME}/fftw-3.3.8/include"
    "$ENV{MAIA_LIB_HOME}/parallel-netcdf-1.8/include"
    "$ENV{MAIA_LIB_HOME}/openmpi/include"
)

# Library directories
set(LIBRARY_DIRS
    "$ENV{MAIA_LIB_HOME}/fftw-3.3.8/lib"
    "$ENV{MAIA_LIB_HOME}/parallel-netcdf-1.8/lib"
    "$ENV{MAIA_LIB_HOME}/openmpi/lib64"
)

# Library names
set(LIBRARY_NAMES
    "fftw3_mpi"
    "fftw3"
    "m"
    "pnetcdf"
    "mpi_cxx"
    "mpi"
)

################################################################################
# Set additional include/library directories or libraries for certain compilers

set(INCLUDE_DIRS_GNU
    "$ENV{MAIA_COMPILER_HOME}/tbb/include" # required by gcc parallel stl
)
# GLIBCXX_3.4.29 required by GCC 11
set(LIBRARY_DIRS_GNU
    "$ENV{MAIA_COMPILER_HOME}/gcc-11.1.0/lib64"
)
set(LIBRARY_NAMES_GNU
    "stdc++"
)

# Additional info for clang's libcxx
set(LIBRARY_DIRS_Clang
    "$ENV{MAIA_COMPILER_HOME}/llvm-9.0.1/lib"
)
set(LIBRARY_NAMES_Clang
    "c++"
    "c++abi"
    "pthread"
)
set(INCLUDE_DIRS_Clang "$ENV{MAIA_COMPILER_HOME}/llvm-9.0.1/include")

# Additional info for Intel's libcxx
set(LIBRARY_DIRS_Intel
    "$ENV{MAIA_COMPILER_HOME}/intel/oneapi/compiler/latest/linux/compiler/include/"
)
set(LIBRARY_NAMES_Intel
    "c++"
    "c++abi"
    "pthread"
)
set(INCLUDE_DIRS_Intel "$ENV{MAIA_COMPILER_HOME}/intel/oneapi/compiler/latest/linux/include")

# NVHPC
set(LIBRARY_DIRS_NVHPC
  "$ENV{MAIA_COMPILER_HOME}/nvhpc-23.7/Linux_x86_64/23.7/compilers/lib"
)

set(INCLUDE_DIRS_NVHPC
  "$ENV{MAIA_COMPILER_HOME}/nvhpc-23.7/Linux_x86_64/23.7/compilers/include"
)

# Support for backtraces with debug(requires Clang)
if(ENABLE_BACKTRACE)
  set(INCLUDE_DIRS ${INCLUDE_DIRS}
      "$ENV{MAIA_COMPILER_HOME}/llvm-9.0.1/include"
  )
  set(LIBRARY_DIRS ${LIBRARY_DIRS}
      "$ENV{MAIA_COMPILER_HOME}/llvm-9.0.1/lib"
  )
  set(LIBRARY_NAMES ${LIBRARY_NAMES}
    "LLVMSupport"
    "LLVMDemangle"
    "dl"
    "tinfo"
    "pthread"
    #"c++"
  )
endif()

#ZLIB
if(WITH_ZLIB)
  set(LIBRARY_NAMES ${LIBRARY_NAMES} "z")
endif()

# Link against HDF5 if enabled
if(WITH_HDF5)
  set(INCLUDE_DIRS ${INCLUDE_DIRS} "$ENV{MAIA_LIB_HOME}/hdf5/include")
  set(LIBRARY_DIRS ${LIBRARY_DIRS} "$ENV{MAIA_LIB_HOME}/hdf5/lib")
  set(LIBRARY_NAMES ${LIBRARY_NAMES} "hdf5_hl")
  set(LIBRARY_NAMES ${LIBRARY_NAMES} "hdf5")
  set(LIBRARY_NAMES ${LIBRARY_NAMES} "z")
  set(LIBRARY_NAMES ${LIBRARY_NAMES} "dl")
endif()

# LIKWID
if(WITH_LIKWID)
  set(INCLUDE_DIRS ${INCLUDE_DIRS} "$ENV{MAIA_LIB_HOME}/likwid/include")
  set(LIBRARY_DIRS ${LIBRARY_DIRS} "$ENV{MAIA_LIB_HOME}/likwid/lib")
  set(LIBRARY_NAMES ${LIBRARY_NAMES} "likwid")
endif()

# CANTERA
if(WITH_CANTERA)
  find_package(Threads REQUIRED)
  set(INCLUDE_DIRS ${INCLUDE_DIRS} ${SRC_DIR_ABS}/../include/canteraSoftLink/include)
  set(LIBRARY_DIRS ${LIBRARY_DIRS} ${SRC_DIR_ABS}/../include/canteraSoftLink/build/lib  )
  set(LIBRARY_NAMES ${LIBRARY_NAMES} "cantera")
  set(CANTERA_DATA ${CANTERA_DATA} ${SRC_DIR_ABS}/../include/canteraSoftLink/build/data)
endif()

#PSTL - cuda version
set(CUDA_VERSION "12.2")

################################################################################
# Enable clang static analyzer for this host
set(CLANG_STATIC_ANALYZER_DIR "$ENV{MAIA_COMPILER_HOME}/llvm-9.0.1/")

################################################################################
# Set (if necessary) configuration commands that should be executed when
# configuring MAIA.
set(HOST_CONFIGURE
    "export LC_ALL=en_US"
    "export LANG=C"
)

set(HOST_CONFIGURE_GNU
    "export LD_LIBRARY_PATH=$ENV{MAIA_COMPILER_HOME}/gcc-11.1.0/lib64:$LD_LIBRARY_PATH")

################################################################################
# Set (if necessary) host-specific default options for the configuration process
# Note: To find out about supported options, check configure.py
set(HOST_DEFAULT_OPTIONS "WITH_HDF5")
set(HOST_DEFAULT_OPTIONS ${HOST_DEFAULT_OPTIONS})
