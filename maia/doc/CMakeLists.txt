cmake_minimum_required(VERSION 2.8.1)

# This script finds all necessary packages and creates a custom target for the
# generation of the Doxygen documentation

################################################################################
# Configure paths - note that if you change any of them, you might need to
# change stuff in the utility scripts or the documentation as well
# Module path for custom CMake modules
set(CMAKE_DIR "cmake")
# Directory with utilities
set(UTILS_DIR "auxiliary")
# Directory with host information files
set(HOSTS_DIR "auxiliary/hosts")
# Set suffix for configuration files for compilers, hosts etc.
set(CONFIG_FILE_SUFFIX ".cmake")
# Directory with solver sources files
set(SOLVER_SRC_DIR "src")
# Directory with documentation files
set(DOC_DIR "doc")

################################################################################
# Set paths
get_filename_component(PARENT_DIR ${PROJECT_SOURCE_DIR} DIRECTORY)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PARENT_DIR}/${CMAKE_DIR})
set(HOSTS_DIR_ABS ${PARENT_DIR}/${HOSTS_DIR})
set(UTILS_DIR_ABS ${PARENT_DIR}/${UTILS_DIR})
set(SOLVER_SRC_DIR_ABS ${PARENT_DIR}/${SOLVER_SRC_DIR})

################################################################################
# Create custom target '${DOC_TARGET}' to generate Doxygen documentation in
# '${DOC_DIR}'
set(_doc "Create 'doc' target to automatically create Doxygen documentation in")
set(_doc "${_doc} '${DOC_DIR}/'.")
set(DOC_TARGET "doc")
set(DOC_TARGET_FAST "doc-fast")
option(ENABLE_DOXYGEN_DOCUMENTATION ${_doc} ON)
unset(_doc)

################################################################################
# Get additional CMake functions
include(GetHost)
include(InList)

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
# Load host file
include(${host_file})

################################################################################
# If documentation is enabled, create documentation target
if(ENABLE_DOXYGEN_DOCUMENTATION)
  # Try to find doxygen
  if(DEFINED HOST_DOXYGEN_EXE)
    set(DOXYGEN_EXECUTABLE ${HOST_DOXYGEN_EXE})
  endif()
  set(Doxygen_FIND_COMPONENTS dot)
  set(_Doxygen_keep_backward_compat TRUE)
  find_package(Doxygen "1.9.5")

  if(DOXYGEN_FOUND)
    # Generate documentation target if doxygen is available
    # Generate project version number if possible
    execute_process(COMMAND
        "git" "describe" "--abbrev=8" "--dirty" "--always" "--tags"
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE project_number
        ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
    if("${project_number}" STREQUAL "")
        set(project_number "n/a")
    endif()

    # Check if DOT was found
    if(DOXYGEN_DOT_FOUND)
      set(have_dot "YES")
    else ()
      set(have_dot "NO")
    endif()

    # Generate Doxygen file
    configure_file(${CMAKE_SOURCE_DIR}/Doxyfile.in
                   ${CMAKE_BINARY_DIR}/Doxyfile @ONLY)

    configure_file(${CMAKE_SOURCE_DIR}/Doxyfile_fast.in
                   ${CMAKE_BINARY_DIR}/Doxyfile_fast @ONLY)

    # Add custom target for the documentation
    add_custom_target(
      ${DOC_TARGET} COMMAND ${DOXYGEN_EXECUTABLE} Doxyfile
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
    add_custom_target(
      ${DOC_TARGET_FAST} COMMAND ${DOXYGEN_EXECUTABLE} Doxyfile_fast
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
  else()
    # Otherwise inform user that doxygen was not found
    add_custom_target(
        ${DOC_TARGET} COMMAND ${UTILS_DIR_ABS}/notfound.py Doxygen)
    add_custom_target(
      ${DOC_TARGET_FAST} COMMAND ${UTILS_DIR_ABS}/notfound.py Doxygen)
  endif()
endif()

