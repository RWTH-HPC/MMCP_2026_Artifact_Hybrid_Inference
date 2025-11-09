################################################################################
# Set available build types
set(BUILD_TYPES "debug" "production" "extreme" "extra_debug" "prototyping"
    "sanitize_address" "sanitize_thread" "sanitize_undefined" "sanitize_conversion")

################################################################################
# Set default build type
set(DEFAULT_BUILD_TYPE "debug")

################################################################################
# Compiler flags

# Default flags
set(_cxx_flags_default "-std=c++17" "-stdlib=libc++" "-Wall" "-Wextra"
    "-pedantic" "-Wshadow" "-Wfloat-equal"
    "-Wcast-align" "-Wfloat-equal" "-Wdisabled-optimization" "-Wformat=2"
    "-Winvalid-pch" "-Winit-self" "-Wmissing-include-dirs" "-Wredundant-decls"
    "-Wpacked" "-Wpointer-arith" "-Wstack-protector" "-Wswitch-default"
    "-Wwrite-strings" "-Wno-type-safety" "-Werror" "-Wunused"
    "-Wno-infinite-recursion" "-fcolor-diagnostics" "-Wno-inconsistent-missing-override"
    # Note: the next warning was added since otherwise too many false-positives
    # are generated. In any case, if a real problem should arise, the linker
    # would find it by complaining about unresolved symbols.
    "-Wno-undefined-var-template" "-Wextra-semi-stmt" )

# Debug
set(_cxx_flags_debug "-O0" "-g3" "-fno-inline" "-DCOMPILER_ATTRIBUTES")

# Production
set(_cxx_flags_production "-O3" "-DNDEBUG" "-fvectorize" "-fslp-vectorize"
    "-DCOMPILER_ATTRIBUTES" "-DUSE_RESTRICT" "-march=native" "-mtune=native")

# Extreme
set(_cxx_flags_extreme "-O3" "-DNDEBUG" "-mtune=native" "-fvectorize"
    "-fslp-vectorize" "-DUSE_RESTRICT" "-DCOMPILER_ATTRIBUTES" "-DDISABLE_FV_MG"
    "-fstrict-aliasing" "-fslp-vectorize-aggressive" "-fno-rtti"
    "-fno-exceptions" "-fomit-frame-pointer" "-march=native" "-mtune=native")

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
set(_cmake_cxx_flags_sanitize_address
    ${_cxx_flags_production}
    "-fsanitize=address" "-DMAIA_SANITIZE_ADDRESS"
    "-g" "-fno-inline" "-fno-omit-frame-pointer"
    "-fno-optimize-sibling-calls"
    ${_cxx_flags_default})
set(_cmake_cxx_flags_sanitize_memory
    ${_cxx_flags_production}
    "-fsanitize=memory" "-fno-inline" "-fno-omit-frame-pointer" "-DMAIA_SANITIZE_MEMORY"
    "-fno-optimize-sibling-calls"
    ${_cxx_flags_default})
set(_cmake_cxx_flags_sanitize_thread
    ${_cxx_flags_production}
    "-fsanitize=thread" "-DMAIA_SANITIZE_THREAD"
    ${_cxx_flags_default})
set(_cmake_cxx_flags_sanitize_undefined
    ${_cxx_flags_production}
    "-fsanitize=undefined,unsigned-integer-overflow,implicit-conversion" "-fno-sanitize=alignment"
    "-fno-sanitize-recover=undefined,unsigned-integer-overflow,implicit-conversion"
    "-DMAIA_SANITIZE_UNDEFINED"
    ${_cxx_flags_default})
set(_cmake_cxx_flags_sanitize_conversion
    ${_cxx_flags_production} "-g3"
    "-fsanitize=implicit-conversion" "-fno-sanitize-recover=implicit-conversion"
    "-DMAIA_SANITIZE_UNDEFINED"
    ${_cxx_flags_default})

join(CMAKE_CXX_FLAGS_DEBUG " " ${_cmake_cxx_flags_debug})
join(CMAKE_CXX_FLAGS_PRODUCTION " " ${_cmake_cxx_flags_production})
join(CMAKE_CXX_FLAGS_EXTREME " " ${_cmake_cxx_flags_extreme})
join(CMAKE_CXX_FLAGS_SANITIZE_ADDRESS " " ${_cmake_cxx_flags_sanitize_address})
join(CMAKE_CXX_FLAGS_SANITIZE_MEMORY " " ${_cmake_cxx_flags_sanitize_memory})
join(CMAKE_CXX_FLAGS_SANITIZE_THREAD " " ${_cmake_cxx_flags_sanitize_thread})
join(CMAKE_CXX_FLAGS_SANITIZE_UNDEFINED " " ${_cmake_cxx_flags_sanitize_undefined})
join(CMAKE_CXX_FLAGS_SANITIZE_CONVERSION " " ${_cmake_cxx_flags_sanitize_conversion})

################################################################################
# Set OpenMP flags (if supported)
set(HAVE_OPENMP "YES")
set(OPENMP_FLAGS "-fopenmp")

################################################################################
# Enable this compiler for the clang static analyzer
set(SUPPORT_CLANG_STATIC_ANALYZER "YES")

################################################################################
# Set flags (if supported) that are used when '-Wunused-xxx' errors should
# be demoted to warnings (if supported)
set(UNUSED_ONLY_WARNING_FLAGS "-Wno-error=unused-parameter"
    "-Wno-error=unused-const-variable" "-Wno-error=unused-result"
    "-Wno-error=unused-function" "-Wno-error=unused-label"
    "-Wno-error=unused-variable" "-Wno-error=unused-value" 
    "-Wno-error=unused-private-field")

################################################################################
# Unset temporary variables
unset(_cxx_flags_default)
unset(_cxx_flags_debug)
unset(_cxx_flags_production)
unset(_cxx_flags_extreme)
unset(_cmake_cxx_flags_debug)
unset(_cmake_cxx_flags_production)
unset(_cmake_cxx_flags_extreme)
unset(_cmake_cxx_flags_sanitize_address)
unset(_cmake_cxx_flags_sanitize_memory)
unset(_cmake_cxx_flags_sanitize_thread)
unset(_cmake_cxx_flags_sanitize_undefined)
unset(_cmake_cxx_flags_sanitize_conversion)
