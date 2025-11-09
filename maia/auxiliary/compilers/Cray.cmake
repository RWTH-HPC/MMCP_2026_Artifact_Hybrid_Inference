################################################################################
# Set available build types
set(BUILD_TYPES "debug" "production" "extreme")

################################################################################
# Set default build type
set(DEFAULT_BUILD_TYPE "production")

################################################################################
# Compiler flags

# Default flags
set(_cxx_flags_default)

# Debug
set(_cxx_flags_debug "-Gn" "-h bounds" "-h dir_check" "-h msglevel_1")

# Production
set(_cxx_flags_production)

# Extreme
set(_cxx_flags_extreme "-O3" "-h display_opt" "-h aggress" "-h cache3"
    "-h ipa5" "-h noomp" "-h thread0" "-h fp3")

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
set(CMAKE_CXX_FLAGS_PRODUCTION "")
join(CMAKE_CXX_FLAGS_EXTREME " " ${_cmake_cxx_flags_extreme})

################################################################################
# Set OpenMP flags (if supported)
set(HAVE_OPENMP "YES")
set(OPENMP_FLAGS "")

################################################################################
# Unset temporary variables
unset(_cxx_flags_default)
unset(_cxx_flags_debug)
unset(_cxx_flags_production)
unset(_cxx_flags_extreme)
unset(_cmake_cxx_flags_debug)
unset(_cmake_cxx_flags_production)
unset(_cmake_cxx_flags_extreme)
