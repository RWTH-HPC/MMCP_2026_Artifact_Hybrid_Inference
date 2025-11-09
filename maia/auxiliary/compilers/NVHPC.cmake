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
#  68: integer conversion resulted in a change of sign
# 111: statement is unreachable
# 128: loop is not reachable (gives false positives)
# 175: subscript out of range
# 177: variable was declared but never referenced
# 185: dynamic initialization in unreachable code
# 497: declaration hides template parameter
# 612: overloaded virtual function only partially overridden
# Feel free to tackle these warnings any time :)
set(_cxx_flags_default "-DCOMPILER_ATTRIBUTES" "--display_error_number"
  "--diag_suppress 68,111,128,175,177,185,497,612"
  "-std=c++17"
)

# Debug
set(_cxx_flags_debug "-O0" "-g")

# Production
#  "-Minfo=all" "-Mneginfo"
set(_cxx_flags_production "-O2" "-DUSE_RESTRICT" "-DNDEBUG")

# Extreme
set(_cxx_flags_extreme "-O4" "-fast" "-Mipa=fast,inline" "-DUSE_RESTRICT"
    "-DNDEBUG")

if(${MAIA_HOST} MATCHES "AIA")
  #TODO: this is a workaround since Eigen library causes problems when using avx512
  list(APPEND _cxx_flags_default "-tp=haswell") # for generic architecture (here, especially very old one)
  list(APPEND _cxx_flags_default "-Mnobuiltin") # to avoid using built-in nvidia math, which is not suited for old architectures
else()
  list(APPEND _cxx_flags_default "-tp=haswell") # # WAR: 3391055
endif()

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
set(OPENMP_FLAGS "-mp")
set(HAVE_PSTL "YES")

################################################################################
# Unset temporary variables
unset(_cxx_flags_default)
unset(_cxx_flags_debug)
unset(_cxx_flags_production)
unset(_cxx_flags_extreme)
unset(_cmake_cxx_flags_debug)
unset(_cmake_cxx_flags_production)
unset(_cmake_cxx_flags_extreme)
