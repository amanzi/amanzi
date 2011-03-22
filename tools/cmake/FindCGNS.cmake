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
#    CGNS_INCLUDE_DIR      (PATH)       Path to the CGNS include file
#    CGNS_INCLUDE_DIRS     (LIST)       List of all required include files
#    CGNS_LIBRARY_DIR      (PATH)       Path to the CGNS library
#    CGNS_LIBRARY          (FILE)       CGNS library
#    CGNS_LIBRARIES        (LIST)       List of all required CGNS libraries
#
#    Additional variables
#    CGNS_VERSION          (STRING)     CGNS Version string
#
# #############################################################################

# Standard CMake modules see CMAKE_ROOT/Modules
include(FindPackageHandleStandardArgs)

# Amanzi CMake functions see <root>/tools/cmake for source
include(PrintVariable)

if ( CGNS_LIBRARIES AND CGNS_INCLUDE_DIRS )

    # Do nothing. Variables are set. No need to search again

else(CGNS_LIBRARIES AND CGNS_INCLUDE_DIRS)

    # Cache variables
    if(CGNS_DIR)
        set(CGNS_DIR "${CGNS_DIR}" CACHE PATH "Path to search for CGNS include and library files")
    endif()

    if(CGNS_INCLUDE_DIR)
        set(CGNS_INCLUDE_DIR "${CGNS_INCLUDE_DIR}" CACHE PATH "Path to search for CGNS include files")
    endif()

    if(CGNS_LIBRARY_DIR)
        set(CGNS_LIBRARY_DIR "${CGNS_LIBRARY_DIR}" CACHE PATH "Path to search for CGNS library files")
    endif()

    
    # Search for include files
    # Search order preference:
    #  (1) CGNS_INCLUDE_DIR - check existence of path AND if the include files exist
    #  (2) CGNS_DIR/<include,include/CGNS++,include/cgns++>
    #  (3) Default CMake paths See cmake --html-help=out.html file for more information.
    #
    set(cgns_inc_names "cgnslib.h")
    if (CGNS_INCLUDE_DIR)

        if (EXISTS "${CGNS_INCLUDE_DIR}")

            find_path(cgns_test_include_path
                      NAMES ${cgns_inc_names}
                      HINTS ${CGNS_INCLUDE_DIR}
                      NO_DEFAULT_PATH)
            if(NOT cgns_test_include_path)
                message(SEND_ERROR "Can not locate ${cgns_inc_names} in ${CGNS_INCLUDE_DIR}")
            endif()
            set(CGNS_INCLUDE_DIR "${cgns_test_include_path}")

        else()
            message(SEND_ERROR "CGNS_INCLUDE_DIR=${CGNS_INCLUDE_DIR} does not exist")
            set(CGNS_INCLUDE_DIR "CGNS_INCLUDE_DIR-NOTFOUND")
        endif()

   else() 

        set(cgns_inc_suffixes "include")
        if(CGNS_DIR)

            if (EXISTS "${CGNS_DIR}" )

                find_path(CGNS_INCLUDE_DIR
                          NAMES ${cgns_inc_names}
                          HINTS ${CGNS_DIR}
                          PATH_SUFFIXES ${cgns_inc_suffixes}
                          NO_DEFAULT_PATH)

            else()
                 message(SEND_ERROR "CGNS_DIR=${CGNS_DIR} does not exist")
                 set(CGNS_INCLUDE_DIR "CGNS_INCLUDE_DIR-NOTFOUND")
            endif()    


        else()

            find_path(CGNS_INCLUDE_DIR
                      NAMES ${cgns_inc_names}
                      PATH_SUFFIXES ${cgns_inc_suffixes})

        endif()

    endif()

    if ( NOT CGNS_INCLUDE_DIR )
        message(SEND_ERROR "Can not locate CGNS include directory")
    endif()

    # Search for libraries 
    # Search order preference:
    #  (1) CGNS_LIBRARY_DIR - check existence of path AND if the include files exist
    #  (2) CGNS_DIR/<lib,lib/CGNS++,lib/cgns++>
    #  (3) Default CMake paths See cmake --html-help=out.html file for more information.
    #
    set(cgns_lib_names "cgns")
    if (CGNS_LIBRARY_DIR)

        if (EXISTS "${CGNS_LIBRARY_DIR}")

            find_library(CGNS_LIBRARY
                         NAMES ${cgns_lib_names}
                         HINTS ${CGNS_LIBRARY_DIR}
                         NO_DEFAULT_PATH)
        else()
            message(SEND_ERROR "CGNS_LIBRARY_DIR=${CGNS_LIBRARY_DIR} does not exist")
            set(CGNS_LIBRARY "CGNS_LIBRARY-NOTFOUND")
        endif()

    else() 

        list(APPEND cgns_lib_suffixes "lib" "Lib")
        if(CGNS_DIR)

            if (EXISTS "${CGNS_DIR}" )

                find_library(CGNS_LIBRARY
                             NAMES ${cgns_lib_names}
                             HINTS ${CGNS_DIR}
                             PATH_SUFFIXES ${cgns_lib_suffixes}
                             NO_DEFAULT_PATH)

            else()
                 message(SEND_ERROR "CGNS_DIR=${CGNS_DIR} does not exist")
                 set(CGNSLIBRARY "CGNS_LIBRARY-NOTFOUND")
            endif()    


        else()

            find_library(CGNS_LIBRARY
                         NAMES ${cgns_lib_names}
                         PATH_SUFFIXES ${cgns_lib_suffixes})

        endif()

    endif()

    if ( NOT CGNS_LIBRARY )
        message(SEND_ERROR "Can not locate CGNS library")
    endif()    

   
    # CGNS does not have any prerequisite libraries
    set(CGNS_INCLUDE_DIRS ${CGNS_INCLUDE_DIR})
    set(CGNS_LIBRARIES    ${CGNS_LIBRARY})

   
endif(CGNS_LIBRARIES AND CGNS_INCLUDE_DIRS )    

# Send useful message if everything is found
find_package_handle_standard_args(CGNS DEFAULT_MSG
                                           CGNS_LIBRARIES
                                           CGNS_INCLUDE_DIRS)

# find_package)handle)standard_args should set CGNS_FOUND but it does not!
if ( CGNS_LIBRARIES AND CGNS_INCLUDE_DIRS)
    set(CGNS_FOUND TRUE)
else()
    set(CGNS_FOUND FALSE)
endif()

# Define the version
if ( CGNS_INCLUDE_DIR )
    set(cgns_h "${CGNS_INCLUDE_DIR}/cgnslib.h")
    file(STRINGS "${cgns_h}" cgns_version_string REGEX "^#define CGNS_VERSION ")
    string(REGEX REPLACE "^#define CGNS_VERSION ([0-9]+).*$" "\\1" cgns_version "${cgns_version_string}")

    #PRINT_VARIABLE(cgns_version_string)
    #PRINT_VARIABLE(cgns_version)

    set(CGNS_VERSION "${cgns_version}")

endif()    

mark_as_advanced(
  CGNS_VERSION
  CGNS_INCLUDE_DIR
  CGNS_INCLUDE_DIRS
  CGNS_LIBRARY
  CGNS_LIBRARIES
  CGNS_LIBRARY_DIR
)
