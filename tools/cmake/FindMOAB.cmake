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
#    MOAB_FOUND            (BOOL)   Flag indicating if MOAB was found
#    MOAB_INCLUDE_DIR      (PATH)   Path to the moab include file
#    MOAB_INCLUDE_DIRS     (LIST)   Paths to include files moab needs
#    MOAB_LIBRARY_DIR      (PATH)   Path to the moab libraries
#    MOAB_LIBRARY          (FILE)   MOAB library
#    MOAB_LIBRARIES        (LIST)   List of libraries to link to MOAB apps
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
                     NAMES MBCore.hpp
                     HINTS ${search_path}
                     PATH_SUFFIXES include
                     DOC "Path to MOAB include files"
                     NO_DEFAULT_PATH)

     if ( NOT MOAB_INCLUDE_DIR ) 
         message(SEND_ERROR "Can not locate MOAB include file in ${MOAB_DIR}")
     endif()

    # Search for library
    find_library( MOAB_LIBRARY
                  NAMES MOAB 
                  HINTS ${search_path}
                  PATH_SUFFIXES lib 
                  DOC "The MOAB library"
                  NO_DEFAULT_PATH)

     if ( MOAB_LIBRARY ) 
         get_filename_component(MOAB_LIBRARY_DIR "${MOAB_LIBRARY}" PATH)
     else()    
         message(SEND_ERROR "Can not locate MOAB library in ${MOAB_DIR}")
     endif()

     # Need logic here to determine if MOAB needs HDF5 and NetCDF
     # For now, I assume (yeah yeah yeah a** u and me) that MOAB
     # requires both. Adding HDF5 libraries to another
     # library string is not a simple append. The FindHDF.cmake
     # module builds HDF5_LIBRARIES with debug and optimized
     # versions. Will need to take this into account when adjusting
     # *_LIBRARIES.
     #find_package(HDF5)
     #set(MOAB_INCLUDE_DIRS ${MOAB_INCLUDE_DIR})
     #list(APPEND MOAB_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS})
     #set(MOAB_LIBRARIES ${MOAB_LIBRARY})
     #list(APPEND MOAB_LIBRARIES ${HDF5_LIBRARIES})

     #find_package(NetCDF)
     set(MOAB_INCLUDE_DIRS ${MOAB_INCLUDE_DIR})
     list(APPEND MOAB_INCLUDE_DIRS ${NetCDF_INCLUDE_DIRS})
     set(MOAB_LIBRARIES ${MOAB_LIBRARY})
     list(APPEND MOAB_LIBRARIES ${NetCDF_LIBRARIES})

     # Remove duplicates
     list(REMOVE_DUPLICATES MOAB_INCLUDE_DIRS)
     list(REMOVE_DUPLICATES MOAB_LIBRARIES)
     PRINT_VARIABLE(MOAB_LIBRARIES)

endif()      

find_package_handle_standard_args(MOAB DEFAULT_MSG
                                  MOAB_LIBRARY
                                  MOAB_INCLUDE_DIR)
mark_as_advanced(
  MOAB_INCLUDE_DIR
  MOAB_INCLUDE_DIRS
  MOAB_LIBRARY
  MOAB_LIBRARY_DIR
  MOAB_LIBRARIES
)

