# -*- mode: cmake -*-

#
# Amanzi CGNS Find Module
#
# Usage:
#    Control the search through CGNS_DIR or setting environment variable
#    CGNS_ROOT to the CGNS installation prefix.
#
#    This module does not search default paths! 
#
#    Following variables are set:
#    CGNS_FOUND            (BOOL)       Flag indicating if CGNS was found
#    CGNS_INCLUDE_DIR      (PATH)       Path to the cgns include files
#    CGNS_LIBRARY_DIR      (PATH)       Path to the cgns libraries
#    CGNS_LIBRARY          (FILE)       CGNS libraries
#
#    Additional variables:
#    CGNS_VERSION          (STRING)     CGNS VERSION STRING
#    
# 
# #############################################################################

# Standard CMake modules see CMAKE_ROOT/Modules
include(FindPackageHandleStandardArgs)

# Amanzi CMake functions see <root>/tools/cmake for source
include(PrintVariable)
include(SelectSearchPath)

# Define the search path 
select_search_path(CGNS search_path search_path_found)
#PRINT_VARIABLE(search_path)
#PRINT_VARIABLE(search_path_found)

if ( search_path_found )
    message(STATUS "Searching for CGNS in ${search_path}")

    set(CGNS_DIR "${search_path}" CACHE PATH "Path to search for CGNS include and library files")

    # Search for include files
    find_path(CGNS_INCLUDE_DIR
                     NAMES cgnslib.h
                     HINTS ${search_path}
                     PATH_SUFFIXES include
                     DOC "Path to CGNS include files"
                     NO_DEFAULT_PATH)

    # Define the version
    if ( CGNS_INCLUDE_DIR )
        set(cgns_h "${CGNS_INCLUDE_DIR}/cgnslib.h")
        file(STRINGS "${cgns_h}" cgns_version_string REGEX "^#define CGNS_VERSION ")
        string(REGEX REPLACE "^#define CGNS_VERSION ([0-9]+).*$" "\\1" cgns_version "${cgns_version_string}")

        #PRINT_VARIABLE(cgns_version_string)
        #PRINT_VARIABLE(cgns_version)

        set(CGNS_VERSION "${cgns_version}")

    endif()    

    # Search for C library
    find_library( CGNS_LIBRARY
                  NAMES cgns
                  HINTS ${search_path}
                  PATH_SUFFIXES lib 
                  DOC "The CGNS library"
                  NO_DEFAULT_PATH)

     if ( NOT CGNS_LIBRARY ) 
         message(SEND_ERROR "Can not locate CGNS library in ${CGNS_DIR}")
     endif()

endif()      

find_package_handle_standard_args(CGNS DEFAULT_MSG
                                  CGNS_LIBRARY
                                  CGNS_INCLUDE_DIR)
mark_as_advanced(
  CGNS_INCLUDE_DIR
  CGNS_C_LIBRARY
  CGNS_VERSION
)

