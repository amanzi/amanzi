# -*- mode: cmake -*-

#
# Amanzi METIS Find Module
#
# Usage:
#    Control the search through METIS_DIR or setting environment variable
#    METIS_ROOT to the METIS installation prefix.
#
#    This module does not search default paths! 
#
#    Following variables are set:
#    METIS_FOUND            (BOOL)       Flag indicating if METIS was found
#    METIS_INCLUDE_DIR      (PATH LIST)  Path to the metis include file
#    METIS_LIBRARY_DIR      (PATH LIST)  Path to the metis libraries
#    METIS_LIBRARY          (FILE)       METIS library
# 
# #############################################################################

# Standard CMake modules see CMAKE_ROOT/Modules
include(FindPackageHandleStandardArgs)

# Amanzi CMake functions see <root>/tools/cmake for source
include(PrintVariable)
include(SelectSearchPath)

# Define the search path 
select_search_path(METIS search_path search_path_found)
#PRINT_VARIABLE(search_path)
#PRINT_VARIABLE(search_path_found)

if ( search_path_found )
    message(STATUS "Searching for METIS in ${search_path}")

    set(METIS_DIR "${search_path}" CACHE PATH "Path to search for METIS include and library files")

    # Search for include files
    find_path(METIS_INCLUDE_DIR
                     NAMES metis.h
                     HINTS ${search_path}
                     PATH_SUFFIXES include
                     DOC "Path to METIS include files"
                     NO_DEFAULT_PATH)

     if ( NOT METIS_INCLUDE_DIR ) 
         message(SEND_ERROR "Can not locate METIS library in ${METIS_DIR}")
     endif()

    # Search for library
    find_library( METIS_LIBRARY
                  NAMES metis
                  HINTS ${search_path}
                  PATH_SUFFIXES lib 
                  DOC "The METIS library"
                  NO_DEFAULT_PATH)

     if ( METIS_LIBRARY ) 
         get_filename_component(METIS_LIBRARY_DIR "${METIS_LIBRARY}" PATH)
     else()    
         message(SEND_ERROR "Can not locate METIS library in ${METIS_DIR}")
     endif()

endif()      

find_package_handle_standard_args(METIS DEFAULT_MSG
                                  METIS_LIBRARY
                                  METIS_INCLUDE_DIR)
mark_as_advanced(
  METIS_INCLUDE_DIR
  METIS_LIBRARY
  METIS_LIBRARY_DIR
)

