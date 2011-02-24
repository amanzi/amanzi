# -*- mode: cmake -*-

#
# Amanzi MOAB Find Module
#
# Usage:
#    Control the search through MOAB_DIR or setting environment variable
#    MOAB_ROOT to the MOAB installation prefix.
#
#    This module does not search default paths! 
#
#    Following variables are set:
#    MOAB_FOUND            (BOOL)       Flag indicating if MOAB was found
#    MOAB_INCLUDE_DIR      (PATH LIST)  Path to the moab include file
#    MOAB_LIBRARY_DIR      (PATH LIST)  Path to the moab libraries
#    MOAB_LIBRARY          (FILE)       MOAB library
# 
# #############################################################################

# Standard CMake modules see CMAKE_ROOT/Modules
include(FindPackageHandleStandardArgs)

# Amanzi CMake functions see <root>/tools/cmake for source
include(PrintVariable)
include(SelectSearchPath)

# Define the search path 
select_search_path(MOAB search_path search_path_found)
#PRINT_VARIABLE(search_path)
#PRINT_VARIABLE(search_path_found)

if ( search_path_found )
    message(STATUS "Searching for MOAB in ${search_path}")

    set(MOAB_DIR "${search_path}" CACHE PATH "Path to search for MOAB include and library files")

    # Search for include files
    find_path(MOAB_INCLUDE_DIR
                     NAMES MBCore.h
                     HINTS ${search_path}
                     PATH_SUFFIXES include
                     DOC "Path to MOAB include files"
                     NO_DEFAULT_PATH)

     if ( NOT MOAB_INCLUDE_DIR ) 
         message(SEND_ERROR "Can not locate MOAB library in ${MOAB_DIR}")
     endif()

    # Search for library
    find_library( MOAB_LIBRARY
                  NAMES MOAB 
                  HINTS ${search_path}
                  PATH_SUFFIXES lib 
                  DOC "The MOAB library"
                  NO_DEFAULT_PATH)

     if ( NOT MOAB_LIBRARY ) 
         message(SEND_ERROR "Can not locate MOAB library in ${MOAB_DIR}")
     endif()

endif()      

find_package_handle_standard_args(MOAB DEFAULT_MSG
                                  MOAB_LIBRARY
                                  MOAB_INCLUDE_DIR)
mark_as_advanced(
  MOAB_INCLUDE_DIR
  MOAB_LIBRARY
)

