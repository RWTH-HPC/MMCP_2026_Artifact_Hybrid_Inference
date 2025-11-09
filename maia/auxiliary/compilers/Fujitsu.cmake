################################################################################
# Set available build types
set(BUILD_TYPES "debug" "production")

################################################################################
# Set default build type
set(DEFAULT_BUILD_TYPE "production")

################################################################################
# Compiler flags

# Default flags
set(_cxx_flags_default "-Xg" "-fPIC" "-std=c++14" "-Bstatic" "-w")

# Debug
set(_cxx_flags_debug "-Ncheck_cache_arraysize"
"-Nquickdbg=subchk" "-Nquickdbg=heapchk" "-Nquickdbg=inf_detail" "-NRtrap")

# Production
set(_cxx_flags_production "-Kfast" )

################################################################################
# Set C++ flags for different build types
set(_cmake_cxx_flags_debug
    ${_cxx_flags_debug}
    ${_cxx_flags_default})
set(_cmake_cxx_flags_production
    ${_cxx_flags_production}
    ${_cxx_flags_default})

join(CMAKE_CXX_FLAGS_DEBUG " " ${_cmake_cxx_flags_debug})
join(CMAKE_CXX_FLAGS_PRODUCTION " " ${_cmake_cxx_flags_production})

################################################################################
# Set OpenMP flags (if supported)
set(HAVE_OPENMP "YES")
set(OPENMP_FLAGS "-Kopenmp")

################################################################################
# Unset temporary variables
unset(_cxx_flags_default)
unset(_cxx_flags_debug)
unset(_cxx_flags_production)
unset(_cmake_cxx_flags_debug)
unset(_cmake_cxx_flags_production)
