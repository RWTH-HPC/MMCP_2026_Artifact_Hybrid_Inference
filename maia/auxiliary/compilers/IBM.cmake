################################################################################
# Set available build types
set(BUILD_TYPES "debug" "production" "extreme")

################################################################################
# Set default build type
set(DEFAULT_BUILD_TYPE "production")

################################################################################
# Compiler flags

# Default flags
set(_cxx_flags_default "-qarch=qp" "-qtune=qp" "-qmaxmem=-1" "-qreport" "-qlist"
    "-qlanglvl=extended0x")

# Debug
set(_cxx_flags_debug "-g" "-qcheck=all" "-qflttrap")

# Production
set(_cxx_flags_production "-O2" "-DNDEBUG")

# Extreme
set(_cxx_flags_extreme "-O4" "-DNDEBUG" "-qstrict")

################################################################################
# Set C++ flags for different build types
set(_cmake_cxx_flags_debug
    ${_cxx_flags_debug}
    ${_cxx_flags_default})
set(_cmake_cxx_flags_production
    ${_cxx_flags_production}
    ${_cxx_flags_default})
set(_cmake_cxx_flags_extreme
    ${_cxx_flags_extreme}
    ${_cxx_flags_default})

join(CMAKE_CXX_FLAGS_DEBUG " " ${_cmake_cxx_flags_debug})
join(CMAKE_CXX_FLAGS_PRODUCTION " " ${_cmake_cxx_flags_production})
join(CMAKE_CXX_FLAGS_EXTREME " " ${_cmake_cxx_flags_extreme})

################################################################################
# Set OpenMP flags (if supported)
set(HAVE_OPENMP "YES")
set(OPENMP_FLAGS "-qsmp")

################################################################################
# Unset temporary variables
unset(_cxx_flags_default)
unset(_cxx_flags_debug)
unset(_cxx_flags_production)
unset(_cxx_flags_extreme)
unset(_cmake_cxx_flags_debug)
unset(_cmake_cxx_flags_production)
unset(_cmake_cxx_flags_extreme)
