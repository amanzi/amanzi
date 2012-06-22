# -*- mode: cmake -*-

#
# Amanzi Mercurial Find Module
#
# Usage:
#
#   Control the search through the variable MERCURIAL_ROOT. Searches for the mercurial
#   binary (i.e. hg) and if it finds the binary defines the version
#
#
#  Module sets
#
#  MERCURIAL_EXECUTABLE   Full path name of the Mercurial binary
#  MERCURIAL_VERSION      Mercurial Version
#  MERCURIAL_FOUND        Flag indicating if mercurial has been found
#
#
include(FindPackageHandleStandardArgs)

if (MERCURIAL_NO_SYSTEM_PATHS)
  set(_hg_FIND_OPTIONS NO_CMAKE_SYSTEM_PATH)
endif()

if ( NOT MERCURIAL_EXECUTABLE )

  find_program(MERCURIAL_EXECUTABLE
               name hg
               PATHS ${MERCURIAL_ROOT}
               ${_hg_FIND_OPTIONS})

endif()

if ( MERCURIAL_EXECUTABLE )

  execute_process(COMMAND ${MERCURIAL_EXECUTABLE} --version
                  RESULT_VARIABLE exit_code
                  OUTPUT_VARIABLE MERCURIAL_VERSION
                  ERROR_VARIABLE  err
                  OUTPUT_STRIP_TRAILING_WHITESPACE
                  ERROR_STRIP_TRAILING_WHITESPACE)

  if (exit_code)
    message(WARNING "Could not determine mercurial version:${err}")
  endif()

endif()

find_package_handle_standard_args(Mercurial DEFAULT_MSG 
                                  MERCURIAL_EXECUTABLE)
                                           


