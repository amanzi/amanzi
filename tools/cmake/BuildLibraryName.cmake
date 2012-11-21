#
# BUILD_LIBRARY_NAME
#
# Usage:
# BUILD_LIBRARY_NAME(library output_name [SHARED|STATIC])
#
# Given a library target name, return the full filename
# of the library with correct suffix and prefix.
# If STATIC or SHARED IS NOT set then the then suffix
# and prefix are defined by BUILD_SHARED_LIBS flag.
# Default is static.
include(CMakeParseArguments)
include(PrintVariable)
function(BUILD_LIBRARY_NAME library output_name)

  set(options "SHARED;STATIC")
  set(oneValue "APPEND_PATH")
  set(multiValue "")
  cmake_parse_arguments(PARSE "${options}" "${oneValue}" "${multiValue}" ${ARGN})

  if ( PARSE_SHARED AND PARSE_STATIC ) 
    message(FATAL_ERROR "Can not ask for STATIC and SHARED library names")
  endif()

  # Set teh suffix and prefix
  set(lib_suffix)
  set(lib_prefix)
  if (PARSE_SHARED OR BUILD_SHARED_LIBS)
    set(lib_suffix ${CMAKE_SHARED_LIBRARY_SUFFIX})
    set(lib_prefix ${CMAKE_SHARED_LIBRARY_PREFIX})
  else ()  
    set(lib_suffix ${CMAKE_STATIC_LIBRARY_SUFFIX})
    set(lib_prefix ${CMAKE_STATIC_LIBRARY_PREFIX})
  endif()

  if ( PARSE_APPEND_PATH )
    set(${output_name} "${PARSE_APPEND_PATH}/${lib_prefix}${library}${lib_suffix}" PARENT_SCOPE)
  else()  
    set(${output_name} "${lib_prefix}${library}${lib_suffix}" PARENT_SCOPE)
  endif()  


endfunction(BUILD_LIBRARY_NAME)
                     
