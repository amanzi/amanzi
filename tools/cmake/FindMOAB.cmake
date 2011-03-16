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
#    MOAB_INCLUDE_DIR      (PATH)       Path to the MOAB include file
#    MOAB_INCLUDE_DIRS     (LIST)       List of all required include files
#    MOAB_LIBRARY_DIR      (PATH)       Path to the MOAB library
#    MOAB_LIBRARY          (FILE)       MOAB library
#    MOAB_LIBRARIES        (LIST)       List of all required MOAB libraries
#
# #############################################################################

# Standard CMake modules see CMAKE_ROOT/Modules
include(FindPackageHandleStandardArgs)

# Amanzi CMake functions see <root>/tools/cmake for source
include(PrintVariable)
include(AddPackageDependency)

if ( MOAB_LIBRARIES AND MOAB_INCLUDE_DIRS )

    # Do nothing. Variables are set. No need to search again

else(MOAB_LIBRARIES AND MOAB_INCLUDE_DIRS)

    # Cache variables
    if(MOAB_DIR)
        set(MOAB_DIR "${MOAB_DIR}" CACHE PATH "Path to search for MOAB include and library files")
    endif()

    if(MOAB_INCLUDE_DIR)
        set(MOAB_INCLUDE_DIR "${MOAB_INCLUDE_DIR}" CACHE PATH "Path to search for MOAB include files")
    endif()

    if(MOAB_LIBRARY_DIR)
        set(MOAB_LIBRARY_DIR "${MOAB_LIBRARY_DIR}" CACHE PATH "Path to search for MOAB library files")
    endif()

    
    # Search for include files
    # Search order preference:
    #  (1) MOAB_INCLUDE_DIR - check existence of path AND if the include files exist
    #  (2) MOAB_DIR/<include>
    #  (3) Default CMake paths See cmake --html-help=out.html file for more information.
    #
    set(moab_inc_names "MBCore.hpp")
    if (MOAB_INCLUDE_DIR)

        if (EXISTS "${MOAB_INCLUDE_DIR}")

            find_path(moab_include_path
                      NAMES ${moab_inc_names}
                      HINTS ${MOAB_INCLUDE_DIR}
                      NO_DEFAULT_PATH)
            if(NOT moab_include_path)
                message(SEND_ERROR "Can not locate ${moab_inc_names} in ${MOAB_INCLUDE_DIR}")
            endif()
            set(MOAB_INCLUDE_DIR "${moab_include_path}")

        else()
            message(SEND_ERROR "MOAB_INCLUDE_DIR=${MOAB_INCLUDE_DIR} does not exist")
            set(MOAB_INCLUDE_DIR "MOAB_INCLUDE_DIR-NOTFOUND")
        endif()

   else() 

        set(moab_inc_suffixes "include")
        if(MOAB_DIR)

            if (EXISTS "${MOAB_DIR}" )

                find_path(MOAB_INCLUDE_DIR
                          NAMES ${moab_inc_names}
                          HINTS ${MOAB_DIR}
                          PATH_SUFFIXES ${moab_inc_suffixes}
                          NO_DEFAULT_PATH)

            else()
                 message(SEND_ERROR "MOAB_DIR=${MOAB_DIR} does not exist")
                 set(MOAB_INCLUDE_DIR "MOAB_INCLUDE_DIR-NOTFOUND")
            endif()    


        else()

            find_path(MOAB_INCLUDE_DIR
                      NAMES ${moab_inc_names}
                      PATH_SUFFIXES ${moab_inc_suffixes})

        endif()

    endif()

    if ( NOT MOAB_INCLUDE_DIR )
        message(SEND_ERROR "Can not locate MOAB include directory")
    endif()

    # Search for libraries 
    # Search order preference:
    #  (1) MOAB_LIBRARY_DIR - check existence of path AND if the library file exists
    #  (2) MOAB_DIR/<lib,Lib>
    #  (3) Default CMake paths See cmake --html-help=out.html file for more information.
    #
    set(moab_lib_names "MOAB")
    if (MOAB_LIBRARY_DIR)

        if (EXISTS "${MOAB_LIBRARY_DIR}")

            find_library(MOAB_LIBRARY
                         NAMES ${moab_lib_names}
                         HINTS ${MOAB_LIBRARY_DIR}
                         NO_DEFAULT_PATH)
        else()
            message(SEND_ERROR "MOAB_LIBRARY_DIR=${MOAB_LIBRARY_DIR} does not exist")
            set(MOAB_LIBRARY "MOAB_LIBRARY-NOTFOUND")
        endif()

    else() 

        list(APPEND moab_lib_suffixes "lib" "Lib")
        if(MOAB_DIR)

            if (EXISTS "${MOAB_DIR}" )

                find_library(MOAB_LIBRARY
                             NAMES ${moab_lib_names}
                             HINTS ${MOAB_DIR}
                             PATH_SUFFIXES ${moab_lib_suffixes}
                             NO_DEFAULT_PATH)

            else()
                 message(SEND_ERROR "MOAB_DIR=${MOAB_DIR} does not exist")
                 set(MOABLIBRARY "MOAB_LIBRARY-NOTFOUND")
            endif()    


        else()

            find_library(MOAB_LIBRARY
                         NAMES ${moab_lib_names}
                         PATH_SUFFIXES ${moab_lib_suffixes})

        endif()

    endif()

    if ( NOT MOAB_LIBRARY )
        message(SEND_ERROR "Can not locate MOAB library")
    endif()    

   
    # Define prerequisite packages
    set(MOAB_INCLUDE_DIRS ${MOAB_INCLUDE_DIR})
    set(MOAB_LIBRARIES    ${MOAB_LIBRARY})
    if (MOAB_NEEDS_NetCDF)
        add_package_dependency(MOAB DEPENDS_ON NetCDF)
    endif()
    if(MOAB_NEEDS_HDF5)
        add_package_dependency(MOAB DEPENDS_ON HDF5)
    endif()    


   
endif(MOAB_LIBRARIES AND MOAB_INCLUDE_DIRS )    

# Send useful message if everything is found
find_package_handle_standard_args(MOAB DEFAULT_MSG
  MOAB_LIBRARIES
  MOAB_INCLUDE_DIRS)

# find_package)handle)standard_args should set MOAB_FOUND but it does not!
if ( MOAB_LIBRARIES AND MOAB_INCLUDE_DIRS)
    set(MOAB_FOUND TRUE)
else()
    set(MOAB_FOUND FALSE)
endif()

# Define the version

mark_as_advanced(
  MOAB_INCLUDE_DIR
  MOAB_INCLUDE_DIRS
  MOAB_LIBRARY
  MOAB_LIBRARIES
  MOAB_LIBRARY_DIR
)
