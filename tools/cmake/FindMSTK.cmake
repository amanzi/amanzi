# -*- mode: cmake -*-

#
# Amanzi MSTK Find Module
#
# Usage:
#    Control the search through MSTK_DIR or setting environment variable
#    MSTK_ROOT to the MSTK installation prefix.
#
#    This module does not search default paths! 
#
#    Following variables are set:
#    MSTK_FOUND            (BOOL)       Flag indicating if MSTK was found
#    MSTK_INCLUDE_DIR      (PATH LIST)  Path to the mstk include file
#    MSTK_LIBRARY_DIR      (PATH LIST)  Path to the mstk libraries
#    MSTK_LIBRARY          (FILE)       MSTK library
# 
# #############################################################################

# Standard CMake modules see CMAKE_ROOT/Modules
include(FindPackageHandleStandardArgs)

# Amanzi CMake functions see <root>/tools/cmake for source
include(PrintVariable)
include(SelectSearchPath)

# Define the search path 
select_search_path(MSTK search_path search_path_found)
#PRINT_VARIABLE(search_path)
#PRINT_VARIABLE(search_path_found)

if ( search_path_found )
    message(STATUS "Searching for MSTK in ${search_path}")

    set(MSTK_DIR "${search_path}" CACHE PATH "Path to search for MSTK include and library files")

    # Search for include files
    find_path(MSTK_INCLUDE_DIR
                     NAMES MSTK.h
                     HINTS ${search_path}
                     PATH_SUFFIXES include
                     DOC "Path to MSTK include files"
                     NO_DEFAULT_PATH)

     if ( NOT MSTK_INCLUDE_DIR ) 
         message(SEND_ERROR "Can not locate MSTK library in ${MSTK_DIR}")
     endif()

    # Search for library
    find_library( MSTK_LIBRARY
                  NAMES mstk
                  HINTS ${search_path}
                  PATH_SUFFIXES lib 
                  DOC "The MSTK library"
                  NO_DEFAULT_PATH)

     if ( MSTK_LIBRARY ) 
         get_filename_component(MSTK_LIBRARY_DIR "${MSTK_LIBRARY}" PATH)
     else()    
         message(SEND_ERROR "Can not locate MSTK library in ${MSTK_DIR}")
     endif()

endif()      

find_package_handle_standard_args(MSTK DEFAULT_MSG
                                  MSTK_LIBRARY
                                  MSTK_INCLUDE_DIR)
mark_as_advanced(
  MSTK_INCLUDE_DIR
  MSTK_LIBRARY
  MSTK_LIBRARY_DIR
)

