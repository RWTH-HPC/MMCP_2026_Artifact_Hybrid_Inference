################################################################################
# Set available build types
set(BUILD_TYPES "debug" "production" "extreme" "sanitize_address"
    "sanitize_memory")

################################################################################
# Set default build type
set(DEFAULT_BUILD_TYPE "production")

################################################################################
# Compiler flags

# Default flags
set(_cxx_flags_default "-std=c++17" 
    "-Wall" "-Wextra" "-Wcomment" "-Wreorder" "-Wreturn-type" "-Wsign-compare" "-Wuninitialized" 
    "-pedantic" "-Wshadow" "-Wfloat-equal"
    "-Wcast-align" "-Wfloat-equal" "-Wformat"
    "-Winvalid-pch" "-Winit-self" "-Wmissing-include-dirs" "-Wredundant-decls"
    "-Wpacked" "-Wpointer-arith" "-Wstack-protector" "-Wswitch-default"
    "-Wwrite-strings" "-Wno-type-safety" "-Werror" "-Wunused-variable" 
    "-Wno-infinite-recursion" "-fcolor-diagnostics" "-Wno-inconsistent-missing-override"
    "-Wno-tautological-constant-compare"
    "-DCOMPILER_ATTRIBUTES"
    # Note: the next warning was added since otherwise too many false-positives
    # are generated. In any case, if a real problem should arise, the linker
    # would find it by complaining about unresolved symbols.
    "-Wno-undefined-var-template"
    "-Wno-non-c-typedef-for-linkage"
    )

# Debug
set(_cxx_flags_debug "-O0" "-g")

# Production
set(_cxx_flags_production "-O2" "-DUSE_RESTRICT" "-DNDEBUG" "-fp-model=precise"
  "-fimf-arch-consistency=true" "-no-fma" "-fimf-arch-consistency=true" )

# Extreme
set(_cxx_flags_extreme "-O3" "-xHOST" "-ipo" "-no-prec-div" "-DUSE_RESTRICT"
    "-DNDEBUG")

# Useful but unsused warning flags (too many errors):
# -Wcheck
# -Wnon-virtual-dtor
# -fast (opt_extreme)

# Add/remove some compiler flags for Julia as they either don't make sense
# (cross compilation!) or because they produce faulty code
if(${MAIA_HOST} MATCHES "Juwels")
  list(REMOVE_ITEM _cxx_flags_extreme "-xHOST")
  list(APPEND _cxx_flags_production "-xMIC-AVX512" "-qopt-report=5" "-qopt-report-phase=vec" "-g")
  list(APPEND _cxx_flags_extreme "-xMIC-AVX512" "-qopt-report=5" "-qopt-report-phase=vec" "-g")  
endif()

# Add/remove some compiler flags for KNL Booster (Jureca) as they either don't make sense
# (cross compilation!) or because they produce faulty code
if(WITH_KNLBOOSTER)
  list(REMOVE_ITEM _cxx_flags_extreme "-xHOST")
  list(APPEND _cxx_flags_production "-xMIC-AVX512")
  list(APPEND _cxx_flags_extreme "-xMIC-AVX512" )  
endif()

if(${MAIA_HOST} MATCHES "Hawk")
  list(APPEND _cxx_flags_default "-march=core-avx2" "-mtune=core-avx2")
endif()

if(${MAIA_HOST} MATCHES "CLAIX")
  list(APPEND _cxx_flags_default "-cxx=icpx")
endif()

if(WITH_CANTERA)
	list(APPEND _cxx_flags_default "-DWITH_CANTERA")	
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
set(_cmake_cxx_flags_sanitize_address
    ${_cxx_flags_debug}
    "-check=stack"
    ${_cxx_flags_default})
set(_cmake_cxx_flags_sanitize_memory
    ${_cxx_flags_debug}
    "-check=stack"
    ${_cxx_flags_default})

join(CMAKE_CXX_FLAGS_DEBUG " " ${_cmake_cxx_flags_debug})
join(CMAKE_CXX_FLAGS_PRODUCTION " " ${_cmake_cxx_flags_production})
join(CMAKE_CXX_FLAGS_EXTREME " " ${_cmake_cxx_flags_extreme})
join(CMAKE_CXX_FLAGS_SANITIZE_ADDRESS " " ${_cmake_cxx_flags_sanitize_address})
join(CMAKE_CXX_FLAGS_SANITIZE_MEMORY " " ${_cmake_cxx_flags_sanitize_memory})

################################################################################
# Set OpenMP flags (if supported)
set(HAVE_OPENMP "YES")
set(OPENMP_FLAGS "-fopenmp")

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
