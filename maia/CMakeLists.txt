cmake_minimum_required(VERSION 2.8.1)

# This is the main CMake configuration file for MAIA

################################################################################
# Set project name, description & set language to C++ (disables C compiler)
project(MAIA CXX)
set(PROJECT_DESCRIPTION "A zonal flow solver." VERSION "")

################################################################################
# Configure paths - note that if you change any of them, you might need to
# change stuff in the utility scripts or the documentation as well
# Module path for custom CMake modules
set(CMAKE_DIR "cmake")
# Output path for binaries
set(BIN_DIR "bin")
# Output path for libraries & archives
set(LIB_DIR "lib")
# Directory with source files
set(SRC_DIR "src")
# Directory with header files
set(INCLUDE_DIR "include")
# Directory with utilities
set(UTILS_DIR "auxiliary")
# Directory with host information files
set(HOSTS_DIR "auxiliary/hosts")
# Directory with compiler information files
set(COMPILERS_DIR "auxiliary/compilers")
# Set suffix for configuration files for compilers, hosts etc.
set(CONFIG_FILE_SUFFIX ".cmake")
# Directory with unit tests
set(TST_DIR "tst")
set(DOCTEST_DIR "include/doctest")

################################################################################
# Configure targets
set(COVERAGE_TARGET "coverage")
set(DISTCLEAN_TARGET "distclean")
set(WHAT_TARGET "what")
set(CLEANSVN_TARGET "clean-svn")
set(COMMIT_TARGET "commit")
set(LIB_TARGET "lib")
set(BIN_TARGET "maia")
set(OBJ_LIB_TARGET "maia_obj_lib")
set(TST_TARGET "test_maia")
#activate to get compiler calls
set(CMAKE_VERBOSE_MAKEFILE OFF)

################################################################################
# Configure output names
set(LIB_NAME "maia")
set(BIN_NAME "maia")


################################################################################
# NO NEED TO CHANGE ANYTHING BEYOND THIS LINE
################################################################################


################################################################################
# Set paths
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/${CMAKE_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${LIB_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${LIB_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${BIN_DIR})
set(SRC_DIR_ABS ${CMAKE_SOURCE_DIR}/${SRC_DIR})
set(INCLUDE_DIR_ABS ${CMAKE_SOURCE_DIR}/${INCLUDE_DIR})
set(TST_DIR_ABS ${CMAKE_SOURCE_DIR}/${TST_DIR})
set(HOSTS_DIR_ABS ${CMAKE_SOURCE_DIR}/${HOSTS_DIR})
set(COMPILERS_DIR_ABS ${CMAKE_SOURCE_DIR}/${COMPILERS_DIR})
set(UTILS_DIR_ABS ${CMAKE_SOURCE_DIR}/${UTILS_DIR})

################################################################################
# Get additional CMake functions
include(GetHost)
include(GetCompiler)
include(GetSubdirectories)
include(InList)
include(Join)

################################################################################
# Check if host is supported
if("$ENV{MAIA_HOST_FILE}" STREQUAL "")
  get_host(MAIA_HOST)
  file(GLOB hosts RELATIVE ${HOSTS_DIR_ABS} ${HOSTS_DIR_ABS}/*)
  in_list(status ${MAIA_HOST}${CONFIG_FILE_SUFFIX} ${hosts})
  if(NOT ${status})
    message(
      FATAL_ERROR "Current host '${MAIA_HOST}' is not recognized. Create a "
      "new file '${MAIA_HOST}${CONFIG_FILE_SUFFIX}' in '${HOSTS_DIR_ABS}' to "
      "support it.")
  endif()
  set(host_file ${HOSTS_DIR_ABS}/${MAIA_HOST}${CONFIG_FILE_SUFFIX})
  set("MAIA_CURRENT_HOST [RO]" "${MAIA_HOST}" CACHE STRING
      "Current host as recognized by MAIA (read-only)." FORCE)
  message(STATUS "Current host as recognized by MAIA: ${MAIA_HOST}.")
else()
  set(host_file $ENV{MAIA_HOST_FILE})
  message(STATUS "Custom host file used: ${host_file}")
endif()

add_definitions(-DHOST_${MAIA_HOST})

################################################################################
# Check if compiler is supported
if(NOT MAIA_COMPILER)
  get_compiler(MAIA_COMPILER)
endif()
file(GLOB compilers RELATIVE ${COMPILERS_DIR_ABS} ${COMPILERS_DIR_ABS}/*)
in_list(status ${MAIA_COMPILER}${CONFIG_FILE_SUFFIX} ${compilers})
if(NOT ${status})
  message(
    FATAL_ERROR
    "Compiler '${MAIA_COMPILER}' not found in '${COMPILERS_DIR_ABS}'."
    " Consider adding a new file '${MAIA_COMPILER}${CONFIG_FILE_SUFFIX}' for it."
    " Alternatively, select a different compiler.")
endif()
set(compiler_file ${COMPILERS_DIR_ABS}/${MAIA_COMPILER}${CONFIG_FILE_SUFFIX})
set("MAIA_COMPILER [RO]" ${MAIA_COMPILER} CACHE
    STRING "Currently selected compiler (read-only)." FORCE)
set("MAIA_COMPILER_VERSION [RO]" ${CMAKE_CXX_COMPILER_VERSION} CACHE
    STRING "Currently selected compiler version (read-only)." FORCE)
message(STATUS "Currently selected compiler: "
        "${MAIA_COMPILER} ${CMAKE_CXX_COMPILER_VERSION}.")

################################################################################
# Load compiler & host files
include(${compiler_file})
include(${host_file})

################################################################################
# Check if compiler is supported on host
join(_supported_compilers ", " ${HOST_SUPPORTED_COMPILERS})
set("MAIA_SUPPORTED_COMPILERS [RO]" ${_supported_compilers} CACHE STRING
    "List of compilers that are supported on the current host (read-only)."
    FORCE
)
unset(_supported_compilers)
in_list(status ${MAIA_COMPILER} ${HOST_SUPPORTED_COMPILERS})
if(NOT ${status})
  message(FATAL_ERROR "Currently selected compiler '${MAIA_COMPILER}' is not "
          "supported on host '${MAIA_HOST}'. Choose a different compiler ("
          "supported: ${_supported_compilers}) or edit '${host_file}'.")
endif()

################################################################################
# Check if selected build type is available
join(_supported_build_types ", " ${BUILD_TYPES})
set("MAIA_SUPPORTED_BUILD_TYPES [RO]" ${_supported_build_types} CACHE STRING
    "List of build types for the selected compiler (read-only)."
    FORCE
)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE ${DEFAULT_BUILD_TYPE})
endif()
if(NOT __INITIALIZED)
set(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING
    "Choose the type of build, options are: ${_supported_build_types}." FORCE)
endif()
message(STATUS "Currently selected build type: ${CMAKE_BUILD_TYPE}")
in_list(status "${CMAKE_BUILD_TYPE}" ${BUILD_TYPES})
if(NOT ${status})
  message(FATAL_ERROR "Build type '${CMAKE_BUILD_TYPE}' is not available for "
          "compiler '${MAIA_COMPILER}' on host '${MAIA_HOST}'. Possible "
          "alternatives are: ${_supported_build_types}")
endif()

################################################################################
# Reset install prefix
if(NOT __INITIALIZED)
  set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR} CACHE PATH
    "Install path prefix, prepended onto install directories." FORCE)
  mark_as_advanced(CMAKE_INSTALL_PREFIX)
endif()

################################################################################
# Add option for creating a MAIA library at all
option(CREATE_LIB "Create MAIA library. Default is `OFF`." OFF)

################################################################################
# Add options for selecting dynamic/static linking
# ...for MAIA itself
option(BUILD_SHARED_LIBS "Build shared MAIA library. Default is `OFF`." OFF)
set(msg "Build static MAIA library, even if `BUILD_SHARED_LIBS` is set to `ON`. "
        "Default is `OFF`.")
join(msg "" ${msg})
option(ALWAYS_BUILD_STATIC ${msg} OFF)
unset(msg)
# ... for libraries linked to MAIA
option(PREFER_STATIC_LIBRARIES
       "Prefer static versions of external linked lobraries." ON)

################################################################################
# Add option for added complete LD_LIBRARY_PATH to the rpath when linking
option(LD_LIBRARY_PATH_TO_RPATH
       "Add all paths in `LD_LIBRARY_PATH' to linker flags." ON)
set(LD_LIBRARY_PATH $ENV{LD_LIBRARY_PATH} CACHE STRING
    "Contents of `LD_LIBRARY_PATH` at configure time. Change only if needed!")

################################################################################
# Process compiler settings
# Add OpenMP flags if OpenMP is enabled
if(ENABLE_OPENMP AND HAVE_OPENMP AND OPENMP_FLAGS)
  string(TOUPPER ${CMAKE_BUILD_TYPE} _bt_upper)
  join(_flags " " ${OPENMP_FLAGS})
  set("CMAKE_CXX_FLAGS_${_bt_upper}"
      "${CMAKE_CXX_FLAGS_${_bt_upper}} ${_flags}")
endif()
# Add parallel STL if enabled
if(ENABLE_PSTL AND HAVE_PSTL AND PSTL_FLAGS)
  string(TOUPPER ${CMAKE_BUILD_TYPE} _bt_upper)
  join(_flags " " ${PSTL_FLAGS})
  set("CMAKE_CXX_FLAGS_${_bt_upper}"
      "${CMAKE_CXX_FLAGS_${_bt_upper}} ${_flags}")
endif()
unset(_bt_upper)
unset(_flags)

# Demote -Wunused-xxx errors to warnings if enabled
if(ALLOW_UNUSED AND UNUSED_ONLY_WARNING_FLAGS)
  string(TOUPPER ${CMAKE_BUILD_TYPE} _bt_upper)
  join(_flags " " ${UNUSED_ONLY_WARNING_FLAGS})
  set("CMAKE_CXX_FLAGS_${_bt_upper}"
      "${CMAKE_CXX_FLAGS_${_bt_upper}} ${_flags}")
endif()
unset(_bt_upper)
unset(_flags)

# Set compiler flags for GUI
foreach(_bt ${BUILD_TYPES})
  string(TOUPPER ${_bt} _bt_upper)
  set(_cmake_cxx_flags_var "CMAKE_CXX_FLAGS_${_bt_upper}")
  mark_as_advanced(${_cmake_cxx_flags_var})
  set(${_cmake_cxx_flags_var} ${${_cmake_cxx_flags_var}} CACHE STRING
      "List of compiler flags for built type '${_bt}'. Change only if needed!")
  if((NOT __INITIALIZED) AND (_bt STREQUAL "debug"))
    set(${_cmake_cxx_flags_var} ${${_cmake_cxx_flags_var}} CACHE STRING
        "List of compiler flags for built type '${_bt}'. Change only if needed!"
        FORCE)
  endif()
endforeach()
unset(_cmake_cxx_flags_var)

################################################################################
# Process host settings
# Include directories
set(INCLUDE_DIRS ${INCLUDE_DIRS} ${INCLUDE_DIRS_${MAIA_COMPILER}})
join(_include_dirs ", " ${INCLUDE_DIRS})
set("MAIA_INCLUDE_DIRS [RO]" "${_include_dirs}" CACHE STRING
    "List of include directories (read-only)."
    FORCE)
mark_as_advanced("MAIA_INCLUDE_DIRS [RO]")
unset(_include_dirs)
# Libraries
if("${PREFER_STATIC_LIBRARIES_PREVIOUS}" STREQUAL "")
  set(PREFER_STATIC_LIBRARIES_PREVIOUS ${PREFER_STATIC_LIBRARIES} CACHE INTERNAL
      "Previous value of static library preference.")
endif()
unset(LIBRARIES)
foreach(_libname ${LIBRARY_NAMES} ${LIBRARY_NAMES_${MAIA_COMPILER}})
  set(_libvar "${_libname}_LIBRARIES")
  if(NOT "${PREFER_STATIC_LIBRARIES}" STREQUAL
     "${PREFER_STATIC_LIBRARIES_PREVIOUS}")
    unset(${_libvar} CACHE)
  endif()
  # If static libraries are preferred, do that only for user-defined libraries,
  # not system libraries
  if(PREFER_STATIC_LIBRARIES)
    set(_lib_searchnames "lib${_libname}.a" ${_libname})
  else()
    set(_lib_searchnames ${_libname})
  endif()
  # Search in user-provided directories first
  find_library(${_libvar} NAMES ${_lib_searchnames} PATHS ${LIBRARY_DIRS}
               ${LIBRARY_DIRS_${MAIA_COMPILER}} NO_DEFAULT_PATH)
  find_library(${_libvar} NAMES ${_libname} PATHS ${LIBRARY_DIRS}
               ${LIBRARY_DIRS_${MAIA_COMPILER}})
  # Search in LD_LIBRARY_PATH (which is for example set by the module environment)
  string(REPLACE ":" ";" RUNTIME_PATH "$ENV{LD_LIBRARY_PATH}")
  find_library(${_libvar} NAMES ${_lib_searchnames} HINTS ${RUNTIME_PATH})
  mark_as_advanced(${_libvar})
  if(${${_libvar}} STREQUAL ${_libvar}-NOTFOUND)
    message(FATAL_ERROR "Could not find library '${_libname}'. "
            "Please check the 'LIBRARY_DIRS' variable in '${host_file}'.")
  endif()
  # Dirty hack but otherwise we get 'ld: cannot find -lc++' errors
  # TODO: The reason for this odd behavior should be found and fixed!
  if(${_libname} STREQUAL "c++")
    get_filename_component(newlib ${${_libvar}} PATH)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${newlib}")
  endif()
  set(LIBRARIES ${LIBRARIES} ${${_libvar}})
endforeach()
set(PREFER_STATIC_LIBRARIES_PREVIOUS ${PREFER_STATIC_LIBRARIES} CACHE INTERNAL
    "Previous value of static library preference." FORCE)

################################################################################
# Create 'analyze' executable when using clang static analyzer
if(CLANG_STATIC_ANALYZER AND CLANG_STATIC_ANALYZER_DIR)
  set(_exe ${HOST_COMPILER_EXE_${MAIA_COMPILER}})
  configure_file(${CMAKE_SOURCE_DIR}/${UTILS_DIR}/analyze.sh.in
                 ${CMAKE_BINARY_DIR}/analyze @ONLY)
  execute_process(COMMAND chmod +x ${CMAKE_BINARY_DIR}/analyze
                  OUTPUT_QUIET ERROR_QUIET)
  unset(_exe)
endif()

################################################################################
# Create custom target '${DISTCLEAN_TARGET}' to remove all files created during
# the CMake build process
add_custom_target(${DISTCLEAN_TARGET}
                  COMMAND ${CMAKE_BUILD_TOOL} clean
                  COMMAND ${CMAKE_SOURCE_DIR}/${UTILS_DIR}/distclean.sh)

################################################################################
# Create custom target to show the current configuration
set(_exe ${HOST_COMPILER_EXE_${MAIA_COMPILER}})
string(TOUPPER ${CMAKE_BUILD_TYPE} _bt_upper)
join(_flags " " ${CMAKE_CXX_FLAGS_${_bt_upper}})
join(_incdirs " " ${INCLUDE_DIRS})
join(_libraries " " ${LIBRARIES})
add_custom_target(
    ${WHAT_TARGET} COMMAND
    ${UTILS_DIR_ABS}/what.py
    ${MAIA_COMPILER}
    ${CMAKE_BUILD_TYPE}
    ${host_file}
    ${compiler_file}
    ${_exe}
    "${_flags}"
    "${_incdirs}"
    "${_libraries}"
    "${LD_LIBRARY_PATH}"
)
unset(_exe)
unset(_bt_upper)
unset(_flags)
unset(_incdirs)
unset(_libraries)

################################################################################
# Create custom targets to clean up as if checked out
add_custom_target(${CLEANSVN_TARGET} COMMAND ${UTILS_DIR_ABS}/cleansvn.sh
                  VERBATIM)
add_custom_target(${CLEANSVN_TARGET}-r COMMAND ${UTILS_DIR_ABS}/cleansvn.sh -r
                  VERBATIM)

################################################################################
# Add targets for preparing a commit
add_custom_target(${COMMIT_TARGET} COMMAND
                  ${UTILS_DIR_ABS}/commit.py --build-system
                  ${CMAKE_SOURCE_DIR}
                  VERBATIM)
add_custom_target(${COMMIT_TARGET}-force COMMAND
                  ${UTILS_DIR_ABS}/commit.py --build-system
                  ${CMAKE_SOURCE_DIR} --force
                  VERBATIM)

################################################################################
# Add targets for analyzing the coverage rwport
add_custom_target(${COVERAGE_TARGET} COMMAND
                  ${UTILS_DIR_ABS}/coverage.py
                  ${SRC_DIR_ABS}
                  ${CMAKE_BINARY_DIR}/${SRC_DIR}/CMakeFiles/maia.dir
                  VERBATIM)
add_custom_target(${COVERAGE_TARGET}-force COMMAND
                  ${UTILS_DIR_ABS}/coverage.py
                  ${SRC_DIR_ABS}
                  ${CMAKE_BINARY_DIR}/${SRC_DIR}/CMakeFiles/maia.dir
                  --force
                  VERBATIM)
add_custom_target(${COVERAGE_TARGET}-verbose COMMAND
                  ${UTILS_DIR_ABS}/coverage.py
                  ${SRC_DIR_ABS}
                  ${CMAKE_BINARY_DIR}/${SRC_DIR}/CMakeFiles/maia.dir
                  -v
                  VERBATIM)

################################################################################
# Add subdirectories with own CMakeLists.txt
add_subdirectory(${SRC_DIR_ABS})
add_subdirectory(${TST_DIR_ABS})
add_subdirectory(${CMAKE_SOURCE_DIR}/${DOCTEST_DIR})

################################################################################
# Set initialization variable - if this is not set, we are during the first
# cmake run and thus may have to override some defaults
set(__INITIALIZED "YES" CACHE INTERNAL "Compiler settings are initialized.")
