################################################################################
# Set available build types
set(BUILD_TYPES "debug" "production" "extreme" "extra_debug" "prototyping"
    "sanitize_address" "sanitize_thread" "coverage" "production_assert"
    "sanitize_undefined")

################################################################################
# Set default build type
set(DEFAULT_BUILD_TYPE "debug")

################################################################################
# Compiler flags

# Default flags
set(_cxx_flags_default "-Werror=return-type" "-Wall" "-Wextra" "-std=c++17" "-pedantic" "-Wshadow"
    "-Wfloat-equal" "-Wcast-align" "-Wfloat-equal" "-Wdisabled-optimization"
    "-Wformat=2" "-Winvalid-pch" "-Winit-self" "-Wmissing-include-dirs"
    "-Wredundant-decls" "-Wpacked" "-Wpointer-arith" "-Wstack-protector"
    "-Wstrict-aliasing=3" "-Wswitch-default" "-Wwrite-strings" "-Wlogical-op"
    "-fdiagnostics-color" "-Wlogical-op" "-Wshift-overflow=2" "-Wnull-dereference"
    "-Wunused-const-variable=1" "-Wduplicated-branches"
    "-Wvla-larger-than=8" "-Walloc-zero"
    "-Wformat-overflow=1" "-Wduplicated-cond" "-Warray-bounds=2"
    "-Wdeprecated-copy-dtor"
    "-Wredundant-tags" "-Wuseless-cast"
    "-Wuninitialized"
    # "-Wmismatched-tags" # gcc-10.2.0 internal compiler error
    #"-Wold-style-cast"
    )

# Removes flag not supported before GCC10
if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 10)
list(REMOVE_ITEM _cxx_flags_default "-Wredundant-tags")
endif()

# Removes flag not supported before GCC9
if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9)
	list(REMOVE_ITEM _cxx_flags_default
		"-Wdeprecated-copy-dtor"
		)
endif()
    
# Removes flag for GCC7 which produces too many false-positives
if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 8)
  list(REMOVE_ITEM _cxx_flags_default
    "-Warray-bounds=2"
    #"-Wold-style-cast"
  )
  list(APPEND _cxx_flags_default
    "-Warray-bounds=0"
  )
endif()

# Remove unsupported flags for older compiler versions
if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 7)
  list(REMOVE_ITEM _cxx_flags_default
    "-Wformat-overflow=1" "-Wduplicated-branches" "-Wduplicated-cond"
    "-Wvla-larger-than=8" "-Walloc-zero" "-Warray-bounds=2"
    #"-Wold-style-cast"
  )
endif()

# Debug
set(_cxx_flags_debug "-Og" "-g3" "-fno-inline" "-gdwarf-4" "-DCOMPILER_ATTRIBUTES"
  "-DMAIA_ASSERTS" "-DMAIA_ASSERT_ALLOC" "-DMAIA_ASSERT_TIMER_CHECKS" "-DMAIA_NCMPI_PRINT_FILE_HINTS"
  "-DMAIA_PROFILING")

# Production
set(_cxx_flags_production "-O3" "-fstrict-aliasing"
    "-fno-rtti" "-fno-exceptions" "-fomit-frame-pointer" 
    "-DCOMPILER_ATTRIBUTES" "-DUSE_RESTRICT" "-DNDEBUG")

set(_cxx_flags_production_assert "-O3" "-g3" "-fstrict-aliasing"
    "-fno-exceptions"
    "-DCOMPILER_ATTRIBUTES" "-DUSE_RESTRICT" "-DMAIA_ASSERTS" "-DMAIA_ASSERT_ACCESSORS")

# Extreme
set(_cxx_flags_extreme ${_cxx_flags_production} "-pipe"
    "-flto=20" "-Wno-unknown-pragmas" "-Wno-maybe-uninitialized"  "-Wno-null-dereference"
    "-fdata-sections" "-ffunction-sections" "-Wl,--gc-sections"
    "-funsafe-loop-optimizations" "-funsafe-math-optimizations"
    "-fcx-limited-range" "-fno-signaling-nans" "-DDISABLE_FV_MG"
    "-fno-math-errno"
    "-Wno-alloc-size-larger-than" # required only for maiamath solveSparseMatrix 
    )

if(${MAIA_HOST} MATCHES "AIA")
  #list(APPEND _cxx_flags_default "-fuse-ld=lld")
  #list(APPEND _cxx_flags_extreme "-fuse-ld=lld") # This is not working !
  list(APPEND _cxx_flags_debug "-fuse-ld=lld")
  list(APPEND _cxx_flags_production "-fuse-ld=lld")
  list(APPEND _cxx_flags_production_assert "-fuse-ld=lld")
endif()

# Add/remove some compiler flags for Power8 as they either don't make sense
# (cross compilation!) or because they produce faulty code
if(${MAIA_HOST} MATCHES "Power8")
  list(REMOVE_ITEM _cxx_flags_production "-march=native" "-mtune=native" "-Wunused")
  list(REMOVE_ITEM _cxx_flags_extreme "-march=native" "-mtune=native" "-Wunused")
endif()

# Remove some compiler flags for GNU compiler on K which are not supported by GCC 4.8.1
if(${MAIA_HOST} MATCHES "K")
  list(REMOVE_ITEM _cxx_flags_default "-fdiagnostics-color" "-Wcast-align")
  list(REMOVE_ITEM _cxx_flags_production "-march=native")
  list(REMOVE_ITEM _cxx_flags_extreme "-march=native")
  list(APPEND _cxx_flags_default "-Wno-missing-field-initializers")
  list(APPEND _cxx_flags_production "-Wno-missing-field-initializers")
  list(APPEND _cxx_flags_extreme "-Wno-missing-field-initializers")
endif()

if(${MAIA_HOST} MATCHES "Juwels")
  list(APPEND _cxx_flags_default "-Wno-useless-cast")
  list(APPEND _cxx_flags_default "-Wno-redundant-tags")
endif()

if(WITH_KNLBOOSTER)
  list(APPEND _cxx_flags_production "-march=knl" "-mtune=knl")
endif()

if(${MAIA_HOST} MATCHES "Hawk")
  list(APPEND _cxx_flags_default "-march=znver2" "-mtune=znver2")
  list(APPEND _cxx_flags_extreme "-Wno-alloc-size-larger-than")
endif()

#CANTERA
if(WITH_CANTERA)
    list(APPEND _cxx_flags_default "-DWITH_CANTERA")
    list(REMOVE_ITEM _cxx_flags_production "-fno-rtti" "-fno-exceptions")
endif()

if(WITH_INSTRUMENTATION)
#<< scorep
# ./configure.py 1 2 --enable-instrumentation --disable-openmp
list(REMOVE_ITEM _cxx_flags_default
  "-Wpedantic"
  "-Wshadow"
  "-Wredundant-decls"
  )
# scorep >>
endif()

################################################################################
# Set C++ flags for different build types
set(_cmake_cxx_flags_debug
    ${_cxx_flags_debug}
    "-Werror"
    ${_cxx_flags_default})
set(_cmake_cxx_flags_production_assert
    ${_cxx_flags_production_assert}
    "-Werror"
    ${_cxx_flags_default})
set(_cmake_cxx_flags_production
    ${_cxx_flags_production}
    "-Werror"
    ${_cxx_flags_default})
set(_cmake_cxx_flags_extreme
    ${_cxx_flags_extreme}
    "-Werror"
    ${_cxx_flags_default})
list(REMOVE_ITEM _cmake_cxx_flags_extreme "-Wnull-dereference")
set(_cmake_cxx_flags_extra_debug
    ${_cxx_flags_debug}
    "-Werror" "-D_GLIBCXX_DEBUG" "-D_FORTIFY_SOURCE=2"
    "-D_GLIBCXX_DEBUG_PEDANTIC" "-Wshift-overflow=2"
    #    "-DMAIA_EXTRA_DEBUG"
    ${_cxx_flags_default})
set(_cmake_cxx_flags_prototyping
    ${_cxx_flags_debug}
    ${_cxx_flags_default})
set(_cmake_cxx_flags_sanitize_address
    ${_cxx_flags_production}
    "-Werror"
    "-fsanitize=address" "-DMAIA_SANITIZE_ADDRESS"
    "-g" "-fno-inline" "-fno-omit-frame-pointer"
    "-fno-optimize-sibling-calls"
    "--param max-gcse-memory=300000000" "--param max-vartrack-size=100000000"
    ${_cxx_flags_default})
set(_cmake_cxx_flags_sanitize_thread
    ${_cxx_flags_production}
    "-Werror"
    "-fsanitize=thread" "-DMAIA_SANITIZE_THREAD"
    ${_cxx_flags_default})
set(_cmake_cxx_flags_sanitize_undefined
    ${_cxx_flags_production}
    "-Werror"
    "-fsanitize=undefined" "-fno-sanitize=alignment" "-DMAIA_SANITIZE_UNDEFINED"
    "-g" "-fno-omit-frame-pointer"
    ${_cxx_flags_default})
set(_cmake_cxx_flags_coverage
    ${_cxx_flags_debug}
    "-Werror"
    "--coverage"
    ${_cxx_flags_default})

# If backtrace is enabled, add flag so that all frames are found
if(ENABLE_BACKTRACE)
  set(_cxx_flags_debug ${_cxx_flags_debug} "-rdynamic")
endif()

# With GCC, the parallel STL requires Intel Thread Building Blocks library
if (ENABLE_PSTL)
  list(APPEND _cxx_flags_default "-ltbb")
endif()

join(CMAKE_CXX_FLAGS_DEBUG " " ${_cmake_cxx_flags_debug})
join(CMAKE_CXX_FLAGS_PRODUCTION_ASSERT " " ${_cmake_cxx_flags_production_assert})
join(CMAKE_CXX_FLAGS_PRODUCTION " " ${_cmake_cxx_flags_production})
join(CMAKE_CXX_FLAGS_EXTREME " " ${_cmake_cxx_flags_extreme})
join(CMAKE_CXX_FLAGS_EXTRA_DEBUG " " ${_cmake_cxx_flags_extra_debug})
join(CMAKE_CXX_FLAGS_PROTOTYPING " " ${_cmake_cxx_flags_prototyping})
join(CMAKE_CXX_FLAGS_SANITIZE_ADDRESS " " ${_cmake_cxx_flags_sanitize_address})
join(CMAKE_CXX_FLAGS_SANITIZE_THREAD " " ${_cmake_cxx_flags_sanitize_thread})
join(CMAKE_CXX_FLAGS_SANITIZE_UNDEFINED " " ${_cmake_cxx_flags_sanitize_undefined})
join(CMAKE_CXX_FLAGS_COVERAGE " " ${_cmake_cxx_flags_coverage})

################################################################################
# Set OpenMP flags (if supported)
set(HAVE_OPENMP "YES")
set(OPENMP_FLAGS "-fopenmp")
set(HAVE_PSTL "YES")

################################################################################
# Enable this compiler for the clang static analyzer
set(SUPPORT_CLANG_STATIC_ANALYZER "YES")

################################################################################
# Set flags (if supported) that are used when '-Wunused-xxx' errors should
# be demoted to warnings (if supported)
set(UNUSED_ONLY_WARNING_FLAGS "-Wno-error=unused-parameter"
    "-Wno-error=unused-but-set-parameter" "-Wno-error=unused-but-set-variable"
    "-Wno-error=unused-function" "-Wno-error=unused-label"
    "-Wno-error=unused-local-typedefs" "-Wno-error=unused-result"
    "-Wno-error=unused-variable" "-Wno-error=unused-value")

################################################################################
# Unset temporary variables
unset(_cxx_flags_default)
unset(_cxx_flags_debug)
unset(_cxx_flags_production_assert)
unset(_cxx_flags_production)
unset(_cxx_flags_extreme)
unset(_cmake_cxx_flags_debug)
unset(_cmake_cxx_flags_production_assert)
unset(_cmake_cxx_flags_production)
unset(_cmake_cxx_flags_extreme)
unset(_cmake_cxx_flags_extra_debug)
unset(_cmake_cxx_flags_prototyping)
unset(_cmake_cxx_flags_sanitize_address)
unset(_cmake_cxx_flags_sanitize_thread)
unset(_cmake_cxx_flags_sanitize_undefined)
unset(_cmake_cxx_flags_coverage)
