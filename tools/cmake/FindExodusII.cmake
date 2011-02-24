# -*- mode: cmake -*-

#
# Amanzi ExodusII Find Module
#
# Usage:
#    Control the search through ExodusII_DIR or setting environment variable
#    ExodusII_ROOT to the ExodusII installation prefix.
#
#    This module does not search default paths! 
#
#    Following variables are set:
#    ExodusII_FOUND            (BOOL)       Flag indicating if ExodusII was found
#    ExodusII_INCLUDE_DIR      (PATH LIST)  Path to the exodus include file
#    ExodusII_LIBRARY_DIR      (PATH LIST)  Path to the exodus libraries
#    ExodusII_LIBRARY          (FILE)       ExodusII library
# 
# #############################################################################

# Standard CMake modules see CMAKE_ROOT/Modules
include(FindPackageHandleStandardArgs)

# Amanzi CMake functions see <root>/tools/cmake for source
include(PrintVariable)
include(SelectSearchPath)

# Define the search path 
select_search_path(ExodusII search_path search_path_found)
#PRINT_VARIABLE(search_path)
#PRINT_VARIABLE(search_path_found)

if ( search_path_found )
    message(STATUS "Searching for ExodusII in ${search_path}")

    set(ExodusII_DIR "${search_path}" CACHE PATH "Path to search for ExodusII include and library files")

    # Search for include files
    find_path(ExodusII_INCLUDE_DIR
                     NAMES exodusII.h
                     HINTS ${search_path}
                     PATH_SUFFIXES include
                     DOC "Path to ExodusII include files"
                     NO_DEFAULT_PATH)

    # Search for library
    find_library( ExodusII_LIBRARY
                  NAMES exoIIv2c 
                  HINTS ${search_path}
                  PATH_SUFFIXES lib 
                  DOC "The ExodusII library"
                  NO_DEFAULT_PATH)

     if ( NOT ExodusII_LIBRARY ) 
         message(SEND_ERROR "Can not locate ExodusII library in ${ExodusII_DIR}")
     endif()

endif()      

find_package_handle_standard_args(ExodusII DEFAULT_MSG
                                  ExodusII_LIBRARY
                                  ExodusII_INCLUDE_DIR)
mark_as_advanced(
  ExodusII_INCLUDE_DIR
  ExodusII_LIBRARY
)

