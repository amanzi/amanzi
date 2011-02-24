# -*- mode: cmake -*-

#
# Amanzi UnitTest Find Module
#
# Usage:
#    Control the search through UnitTest_DIR or setting environment variable
#    UnitTest_ROOT to the UnitTest installation prefix.
#
#    This module does not search default paths! 
#
#    Following variables are set:
#    UnitTest_FOUND            (BOOL)       Flag indicating if UnitTest was found
#    UnitTest_INCLUDE_DIR      (PATH)       Path to the utest include files
#    UnitTest_LIBRARY_DIR      (PATH)       Path to the utest libraries
#    UnitTest_LIBRARY          (FILE)       UnitTest libraries
#
# #############################################################################

# Standard CMake modules see CMAKE_ROOT/Modules
include(FindPackageHandleStandardArgs)

# Amanzi CMake functions see <root>/tools/cmake for source
include(PrintVariable)
include(SelectSearchPath)

# Define the search path 
select_search_path(UnitTest search_path search_path_found)
#PRINT_VARIABLE(search_path)
#PRINT_VARIABLE(search_path_found)

if ( search_path_found )
    message(STATUS "Searching for UnitTest in ${search_path}")

    set(UnitTest_DIR "${search_path}" CACHE PATH "Path to search for UnitTest include and library files")

    # Search for include files
    find_path(UnitTest_INCLUDE_DIR
                     NAMES UnitTest++.h unittest++.h
                     HINTS ${search_path}
                     PATH_SUFFIXES include
                     DOC "Path to UnitTest include files"
                     NO_DEFAULT_PATH)

    # Search for C library
    find_library( UnitTest_LIBRARY
                  NAMES UnitTest++ unittest++
                  HINTS ${search_path}
                  PATH_SUFFIXES lib 
                  DOC "The UnitTest library"
                  NO_DEFAULT_PATH)

     if ( NOT UnitTest_LIBRARY ) 
         message(SEND_ERROR "Can not locate UnitTest library in ${UnitTest_DIR}")
     endif()

endif()      

find_package_handle_standard_args(UnitTest DEFAULT_MSG
                                  UnitTest_LIBRARY
                                  UnitTest_INCLUDE_DIR)
mark_as_advanced(
  UnitTest_INCLUDE_DIR
  UnitTest_LIBRARY
  UnitTest_VERSION
)

