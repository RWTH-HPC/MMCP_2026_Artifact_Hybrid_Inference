################################################################################
# Set available build types
set(BUILD_TYPES "debug" "production" "extreme")

################################################################################
# Set default build type
set(DEFAULT_BUILD_TYPE "production")

################################################################################
# Compiler flags

# Default flags
# The following warnings are suppressed:
#   1: last line of file ends without a newline
# 128: loop is not reachable (gives false positives)
# 185: dynamic initialization in unreachable code
# 490: template cannot be instantiated -- it has been explicitly specialized
# 611: overloaded virtual function only partially overridden
# Feel free to tackle these warnings any time :)
set(_cxx_flags_default "-DCOMPILER_ATTRIBUTES" "--display_error_number"
    "--diag_suppress 1,128,185,490,611" "--c++17")

# Debug
set(_cxx_flags_debug "-O0" "-g" "-Minfo=all" "-Mneginfo")

# Production
set(_cxx_flags_production "-O2" "-DUSE_RESTRICT" "-DNDEBUG")

# Extreme
set(_cxx_flags_extreme "-O4" "-fast" "-Mipa=fast,inline" "-DUSE_RESTRICT"
    "-DNDEBUG")

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
set(HAVE_OPENMP "NO")
set(OPENMP_FLAGS "-mp")

################################################################################
# Unset temporary variables
unset(_cxx_flags_default)
unset(_cxx_flags_debug)
unset(_cxx_flags_production)
unset(_cxx_flags_extreme)
unset(_cmake_cxx_flags_debug)
unset(_cmake_cxx_flags_production)
unset(_cmake_cxx_flags_extreme)
