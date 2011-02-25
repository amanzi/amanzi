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

     if ( ExodusII_LIBRARY ) 
         get_filename_component(ExodusII_LIBRARY_DIR "${ExodusII_LIBRARY}" PATH)
     else()    
         message(SEND_ERROR "Can not locate ExodusII library in ${ExodusII_DIR}")
     endif()

     # ExodusII requires NetCDF. Add extra files and paths to ExodusII_INCLUDE_DIRS
     # and ExodusII_LIBRARIES
     set(ExodusII_INCLUDE_DIRS ${ExodusII_INCLUDE_DIR})
     set(ExodusII_LIBRARIES ${ExodusII_LIBRARY})
     if ( NetCDF_INCLUDE_DIRS AND NetCDF_LIBRARIES )
        # Do nothing
     else ()
         find_package(NetCDF REQUIRED)
     endif()
     list(APPEND ExodusII_INCLUDE_DIRS "${NetCDF_INCLUDE_DIRS}")
     list(APPEND ExodusII_LIBRARIES    "${NetCDF_LIBRARIES}")
     #if ( ${NetCDF_FOUND} )
     #    list(APPEND ExodusII_INCLUDE_DIRS "${NetCDF_INCLUDE_DIRS}")
     #    list(APPEND ExodusII_LIBRARIES    "${NetCDF_LIBRARIES}")
     #else()
     #    message(STATUS "Can not locate NetCDF! This is a required package for ExodusII")
     #endif()

     # Remove duplicates
     list(REMOVE_DUPLICATES ExodusII_INCLUDE_DIRS)
     list(REMOVE_DUPLICATES ExodusII_LIBRARIES)

endif()      

find_package_handle_standard_args(ExodusII DEFAULT_MSG
                                  ExodusII_LIBRARY
                                  ExodusII_INCLUDE_DIR)
mark_as_advanced(
  ExodusII_INCLUDE_DIR
  ExodusII_LIBRARY
)

