# -*- mode: cmake -*-

#
# Amanzi NetCDF Find Module
#
# Usage:
#    Control the search through NetCDF_DIR or setting environment variable
#    NetCDF_ROOT to the NetCDF installation prefix.
#
#    This module does not search default paths! 
#
#    Following variables are set:
#    NetCDF_FOUND            (BOOL)       Flag indicating if NetCDF was found
#    NetCDF_INCLUDE_DIRS     (PATH LIST)  Path to the netcdf include files
#    NetCDF_INCLUDE_DIR      (PATH LIST)  Path to the netcdf include files (deprecated)
#    NetCDF_LIBRARY_DIRS     (PATH LIST)  Path to the netcdf libraries
#    NetCDF_LIBRARIES        (FILE LIST)  NetCDF libraries
# 
#    Additional variables set
#    NetCDF_C_LIBRARY        (FILE)       NetCDF C library
#    NetCDF_CXX_LIBRARY      (FILE)       NetCDF C++ library
#    NetCDF_FORTRAN_LIBRARY  (FILE)       NetCDF Fortran library (only set if Fortran enabled)
#    NetCDF_LARGE_DIMS       (BOOL)       Checks the header files for size of 
#                                          NC_MAX_DIMS, NC_MAX_VARS and NC_MAX_VARS_DIMS
#                                          Returns TRUE if
#                                          NC_MAX_DIMS >= 655363
#                                          NC_MAX_VARS >= 524288
#                                          NC_MAX_VAR_DIMS >= 8
#
# #############################################################################

# Standard CMake modules see CMAKE_ROOT/Modules
include(FindPackageHandleStandardArgs)

# Amanzi CMake functions see <root>/tools/cmake for source
include(PrintVariable)
include(SelectSearchPath)

# Define the search path 
select_search_path(NetCDF search_path search_path_found)
#PRINT_VARIABLE(search_path)
#PRINT_VARIABLE(search_path_found)

if ( search_path_found )
    message(STATUS "Searching for NetCDF in ${search_path}")

    set(NetCDF_DIR "${search_path}" CACHE PATH "Path to search for NetCDF include and library files")

    # Search for include files
    find_path(NetCDF_INCLUDE_DIRS
                     NAMES netcdf.h
                     HINTS ${search_path}
                     PATH_SUFFIXES include
                     DOC "Path to NetCDF include files"
                     NO_DEFAULT_PATH)
    set(NetCDF_INCLUDE_DIR "${NetCDF_INCLUDE_DIRS}")  # Consistent with other standard CMake Modules             

    # MOAB check here
    if ( NetCDF_INCLUDE_DIRS ) 
       
        set(netcdf_h "${NetCDF_INCLUDE_DIRS}/netcdf.h" )

        file(STRINGS "${netcdf_h}" netcdf_max_dims_string REGEX "^#define NC_MAX_DIMS ")
        string(REGEX REPLACE "^#define NC_MAX_DIMS ([0-9]+).*$" "\\1" netcdf_max_dims "${netcdf_max_dims_string}")

        file(STRINGS "${netcdf_h}" netcdf_max_vars_string REGEX "^#define NC_MAX_VARS ")
        string(REGEX REPLACE "^#define NC_MAX_VARS ([0-9]+).*$" "\\1" netcdf_max_vars "${netcdf_max_vars_string}")

        file(STRINGS "${netcdf_h}" netcdf_max_var_dims_string REGEX "^#define NC_MAX_VAR_DIMS ")
        string(REGEX REPLACE "^#define NC_MAX_VAR_DIMS ([0-9]+).*$" "\\1" netcdf_max_var_dims "${netcdf_max_var_dims_string}")

        #PRINT_VARIABLE(netcdf_max_dims)
        #PRINT_VARIABLE(netcdf_max_vars)
        #PRINT_VARIABLE(netcdf_max_var_dims)

        if ( ( netcdf_max_dims EQUAL "65536" OR netcdf_max_dims GREATER "65536") AND
             ( netcdf_max_vars EQUAL "524288" OR netcdf_max_vars GREATER "524288") AND
             ( netcdf_max_var_dims EQUAL "8" OR netcdf_max_var_dims GREATER "8")
             )
            set(NetCDF_MOAB_COMPATIBLE TRUE)
        else()
            message(WARNING "The NetCDF found in ${NetCDF_DIR} does not have the correct NC_MAX_DIMS, NC_MAX_VARS and NC_MAX_VAR_DIMS\n"
                             "It may not be compatible with other TPL libraries such MOAB and ExodusII\n" )
            set(NetCDF_MOAB_COMPATIBLE FALSE)
        endif()

    endif()    

    # Search for C library
    find_library( NetCDF_C_LIBRARY
                  NAMES netcdf
                  HINTS ${search_path}
                  PATH_SUFFIXES lib 
                  DOC "The NetCDF C library"
                  NO_DEFAULT_PATH)

     if ( NOT NetCDF_C_LIBRARY ) 
         message(SEND_ERROR "Can not locate NetCDF C library in ${NetCDF_DIR}")
     endif()

     #Search for C++ library
     find_library(NetCDF_CXX_LIBRARY netcdf_c++
                  HINTS ${search_path}
                  PATH_SUFFIXES lib
                  DOC "The NetCDF C++ library"
                  NO_DEFAULT_PATH)

     if ( NOT NetCDF_CXX_LIBRARY ) 
         message(SEND_ERROR "Can not locate NetCDF C++ library in ${NetCDF_DIR}")
     endif()

     set(NetCDF_LIBRARIES
                         ${NetCDF_C_LIBRARY}
                         ${NetCDF_CXX_LIBRARY} 
                         )

     # Search Fortran library
     if ( ENABLE_FORTRAN )
         find_library(NetCDF_FORTRAN_LIBRARY
                      NAMES netcdf_g77 netcdf_ifc netcdf_x86_64
                      HINTS ${search_path}
                      PATH_SUFFIXES lib
                      DOC "The NetCDF C++ library"
                      NO_DEFAULT_PATH)

         if ( NOT NetCDF_FORTRAN_LIBRARY ) 
             message(SEND_ERROR "Can not locate NetCDF Fortran library in ${NetCDF_DIR}")
         endif()

      endif()

endif()      

find_package_handle_standard_args(NetCDF DEFAULT_MSG
                                  NetCDF_LIBRARIES
                                  NetCDF_INCLUDE_DIRS)
mark_as_advanced(
  NetCDF_INCLUDE_DIR
  NetCDF_C_LIBRARY
  NetCDF_CXX_LIBRARY
  NetCDF_FORTRAN_LIBRARY
  NetCDF_MOAB_COMPATIBLE
)

