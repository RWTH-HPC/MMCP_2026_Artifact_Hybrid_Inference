################################################################################
# Include helper functions
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${MAIA_ROOT}/${CMAKE_DIR})
include(Join)
include(GetHost)
include(InList)

################################################################################
# Add own helpers
function(echo _msg)
  message(${_msg})
endfunction()

################################################################################
# Set paths
set(HOSTS_DIR_ABS ${MAIA_ROOT}/${HOSTS_DIR})
set(COMPILERS_DIR_ABS ${MAIA_ROOT}/${COMPILERS_DIR})

################################################################################
# Get list of available hosts
file(GLOB hosts RELATIVE
     ${HOSTS_DIR_ABS} ${HOSTS_DIR_ABS}/*${CONFIG_FILE_SUFFIX})

################################################################################
# Get and return host
get_host(MAIA_HOST)
if("$ENV{MAIA_HOST_FILE}" STREQUAL "")
  in_list(status ${MAIA_HOST}${CONFIG_FILE_SUFFIX} ${hosts})
  if(status)
    echo("HOST ${MAIA_HOST}")
  else()
    echo("HOST ${MAIA_HOST} NOTFOUND")
    return()
  endif()
  set(host_file ${HOSTS_DIR_ABS}/${MAIA_HOST}${CONFIG_FILE_SUFFIX})
else()
  echo("HOST ${MAIA_HOST}")
  set(host_file $ENV{MAIA_HOST_FILE})
endif()

################################################################################
# Get and print list of available compilers on this host
set(MAIA_COMPILER "UNKNOWN")
include(${host_file})
set(compilers ${HOST_SUPPORTED_COMPILERS})
join(compilersstring " " ${compilers})
echo("COMPILERS ${compilersstring}")
echo("DEFAULT_COMPILER ${HOST_DEFAULT_COMPILER}")
foreach(_compiler ${compilers})
  set(varname "HOST_COMPILER_EXE_${_compiler}")
  echo("COMPILER_EXE ${_compiler} ${${varname}}")
endforeach()

# Get and return required cuda version if available
if(NOT "${CUDA_VERSION}" STREQUAL "")
  echo("CUDA_VERSION ${CUDA_VERSION}")
endif()

################################################################################
# Print info for clang static analyzer if available
if(CLANG_STATIC_ANALYZER_DIR)
  echo("CLANG_STATIC_ANALYZER_DIR ${CLANG_STATIC_ANALYZER_DIR}")
endif()

################################################################################
# Get and return supported optimization levels for each compiler
foreach(compiler ${compilers})
  set(compiler_file ${COMPILERS_DIR_ABS}/${compiler}${CONFIG_FILE_SUFFIX})
  include(${compiler_file})
  join(build_types " " ${BUILD_TYPES})
  echo("BUILD_TYPES ${compiler} ${build_types}")
  echo("DEFAULT_BUILD_TYPE ${compiler} ${DEFAULT_BUILD_TYPE}")
  if(HAVE_OPENMP)
    echo("HAVE_OPENMP ${compiler}")
    join(openmp_flags " " ${OPENMP_FLAGS})
    echo("OPENMP_FLAGS ${compiler} ${openmp_flags}")
    unset(HAVE_OPENMP)
    unset(OPENMP_FLAGS)
  endif()
  if(HAVE_PSTL)
    echo("HAVE_PSTL ${compiler}")
    join(pstl_flags " " ${PSTL_FLAGS})
    echo("PSTL_FLAGS ${compiler} ${pstl_flags}")
    unset(HAVE_PSTL)
    unset(PSTL_FLAGS)
  endif()
  if(SUPPORT_CLANG_STATIC_ANALYZER)
    echo("SUPPORT_CLANG_STATIC_ANALYZER ${compiler}")
    unset(SUPPORT_CLANG_STATIC_ANALYZER)
  endif()
endforeach()

################################################################################
# Return configure commands for each compiler
if(DEFINED HOST_CONFIGURE)
  foreach(_c ${HOST_CONFIGURE})
    echo("CONFIGURE_HOST ${_c}")
  endforeach()
endif()
foreach(compiler ${compilers})
  set(varname "HOST_CONFIGURE_${compiler}")
  if(DEFINED ${varname})
    foreach(_c ${${varname}})
      echo("CONFIGURE ${compiler} ${_c}")
  endforeach()
  endif()
endforeach()

################################################################################
# Return environment variables that should be preserved
if(DEFINED HOST_PRESERVE_ENV)
  join(_vars " " ${HOST_PRESERVE_ENV})
  echo("PRESERVE_ENV_HOST ${_vars}")
endif()
foreach(compiler ${compilers})
  set(varname "HOST_PRESERVE_ENV_${compiler}")
  if(DEFINED ${varname})
    foreach(_c ${${varname}})
      echo("PRESERVE_ENV ${compiler} ${_c}")
  endforeach()
  endif()
endforeach()

################################################################################
# Return generic default options
if(DEFINED HOST_DEFAULT_OPTIONS)
  join(_vars " " ${HOST_DEFAULT_OPTIONS})
  echo("DEFAULT_OPTIONS ${_vars}")
endif()
