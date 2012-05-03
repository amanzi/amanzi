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
#    ExodusII_INCLUDE_DIR      (PATH)       Path to the ExodusII include file
#    ExodusII_INCLUDE_DIRS     (LIST)       List of all required include files
#    ExodusII_LIBRARY_DIR      (PATH)       Path to the ExodusII library
#    ExodusII_LIBRARY          (FILE)       ExodusII library
#    ExodusII_LIBRARIES        (LIST)       List of all required ExodusII libraries
#
#    Additional variables
#    ExodusII_VERSION          (STRING)     ExodusII Version string
#
# #############################################################################

# Standard CMake modules see CMAKE_ROOT/Modules
include(FindPackageHandleStandardArgs)

include(AddImportedLibrary)

# Amanzi CMake functions see <root>/tools/cmake for source
include(PrintVariable)
include(AddPackageDependency)

if ( ExodusII_LIBRARIES AND ExodusII_INCLUDE_DIRS )

    # Do nothing. Variables are set. No need to search again

else(ExodusII_LIBRARIES AND ExodusII_INCLUDE_DIRS)

    # Cache variables
    if(ExodusII_DIR)
        set(ExodusII_DIR "${ExodusII_DIR}" CACHE PATH "Path to search for ExodusII include and library files")
    endif()

    if(ExodusII_INCLUDE_DIR)
        set(ExodusII_INCLUDE_DIR "${ExodusII_INCLUDE_DIR}" CACHE PATH "Path to search for ExodusII include files")
    endif()

    if(ExodusII_LIBRARY_DIR)
        set(ExodusII_LIBRARY_DIR "${ExodusII_LIBRARY_DIR}" CACHE PATH "Path to search for ExodusII library files")
    endif()

    
    # Search for include files
    # Search order preference:
    #  (1) ExodusII_INCLUDE_DIR - check existence of path AND if the include files exist
    #  (2) ExodusII_DIR/<include>
    #  (3) Default CMake paths See cmake --html-help=out.html file for more information.
    #
    set(exodus_inc_names "exodusII.h")
    if (ExodusII_INCLUDE_DIR)

        if (EXISTS "${ExodusII_INCLUDE_DIR}")

            find_path(exodusII_test_include_path
                      NAMES ${exodus_inc_names}
                      HINTS ${ExodusII_INCLUDE_DIR}
                      NO_DEFAULT_PATH)
            if(NOT exodusII_test_include_path)
                message(SEND_ERROR "Can not locate ${exodus_inc_names} in ${ExodusII_INCLUDE_DIR}")
            endif()
            set(ExodusII_INCLUDE_DIR "${exodusII_test_include_path}")

        else()
            message(SEND_ERROR "ExodusII_INCLUDE_DIR=${ExodusII_INCLUDE_DIR} does not exist")
            set(ExodusII_INCLUDE_DIR "ExodusII_INCLUDE_DIR-NOTFOUND")
        endif()

   else() 

        set(exodus_inc_suffixes "include")
        if(ExodusII_DIR)

            if (EXISTS "${ExodusII_DIR}" )

                find_path(ExodusII_INCLUDE_DIR
                          NAMES ${exodus_inc_names}
                          HINTS ${ExodusII_DIR}
                          PATH_SUFFIXES ${exodus_inc_suffixes}
                          NO_DEFAULT_PATH)

            else()
                 message(SEND_ERROR "ExodusII_DIR=${ExodusII_DIR} does not exist")
                 set(ExodusII_INCLUDE_DIR "ExodusII_INCLUDE_DIR-NOTFOUND")
            endif()    


        else()

            find_path(ExodusII_INCLUDE_DIR
                      NAMES ${exodus_inc_names}
                      PATH_SUFFIXES ${exodus_inc_suffixes})

        endif()

    endif()

    if ( NOT ExodusII_INCLUDE_DIR )
        message(SEND_ERROR "Can not locate ExodusII include directory")
    endif()

    # Search for libraries 
    # Search order preference:
    #  (1) ExodusII_LIBRARY_DIR - check existence of path AND if the library file exists
    #  (2) ExodusII_DIR/<lib,Lib>
    #  (3) Default CMake paths See cmake --html-help=out.html file for more information.
    #
    set(exodus_lib_names "exoIIv2c" "exodus")
    if (ExodusII_LIBRARY_DIR)

        if (EXISTS "${ExodusII_LIBRARY_DIR}")

            find_library(_ExodusII_LIBRARY
                         NAMES ${exodus_lib_names}
                         HINTS ${ExodusII_LIBRARY_DIR}
                         NO_DEFAULT_PATH)

            find_library(_ExodusII_Fortran_LIBRARY
                         NAMES exodus_for
                         HINTS ${ExodusII_LIBRARY_DIR}
                         NO_DEFAULT_PATH)

        else()
            message(SEND_ERROR "ExodusII_LIBRARY_DIR=${ExodusII_LIBRARY_DIR} does not exist")
            set(_ExodusII_LIBRARY "ExodusII_LIBRARY-NOTFOUND")
            set(_ExodusII_Fortran_LIBRARY "ExodusII_Fortran_LIBRARY-NOTFOUND")
        endif()

    else() 

        list(APPEND exodus_lib_suffixes "lib" "Lib")
        if(ExodusII_DIR)

            if (EXISTS "${ExodusII_DIR}" )

                find_library(_ExodusII_LIBRARY
                             NAMES ${exodus_lib_names}
                             HINTS ${ExodusII_DIR}
                             PATH_SUFFIXES ${exodus_lib_suffixes}
                             NO_DEFAULT_PATH)
                
		find_library(_ExodusII_Fortran_LIBRARY
                             NAMES exodus_for
                             HINTS ${ExodusII_DIR}
                             PATH_SUFFIXES ${exodus_lib_suffixes}
                             NO_DEFAULT_PATH)

            else()
                 message(SEND_ERROR "ExodusII_DIR=${ExodusII_DIR} does not exist")
                 set(ExodusII_LIBRARY "ExodusII_LIBRARY-NOTFOUND")
                 set(ExodusII_Fortran_LIBRARY "ExodusII_Fortran_LIBRARY-NOTFOUND")
            endif()    


        else()

            find_library(_ExodusII_LIBRARY
                         NAMES ${exodus_lib_names}
                         PATH_SUFFIXES ${exodus_lib_suffixes})

            find_library(_ExodusII_Fortran_LIBRARY
                         NAMES exodus_for
                         PATH_SUFFIXES ${exodus_lib_suffixes})

        endif()

    endif()

    # Create the library target store the name in ExodusII_LIBRARY
    if ( _ExodusII_LIBRARY )
        set(ExodusII_LIBRARY exodusii)
        add_imported_library(${ExodusII_LIBRARY}
                     LOCATION ${_ExodusII_LIBRARY}
                     LINK_LANGUAGES "C;CXX")
    else()
        message(SEND_ERROR "Can not locate ExodusII library")
    endif()

    if ( _ExodusII_Fortran_LIBRARY )
        set(ExodusII_Fortran_LIBRARY exodusii_for)
        add_imported_library(${ExodusII_Fortran_LIBRARY}
                             LOCATION ${_ExodusII_Fortran_LIBRARY}
                             LINK_LANGUAGES "Fortran")
        set_target_properties(${ExodusII_Fortran_LIBRARY} PROPERTIES
	                      IMPORTED_LINK_INTERFACE_LIBRARIES "${Exodus_LIBRARY}")
    endif()


   
    # Define prerequisite packages
    set(ExodusII_INCLUDE_DIRS ${ExodusII_INCLUDE_DIR})
    set(ExodusII_LIBRARIES    ${ExodusII_LIBRARY})

    # Search for NetCDF
    find_package(NetCDF QUIET REQUIRED)
    set_target_properties(${ExodusII_LIBRARY} PROPERTIES
                          IMPORTED_LINK_INTERFACE_LIBRARIES "${NetCDF_LIBRARIES}")
    list(APPEND ExodusII_INCLUDE_DIRS ${NetCDF_INCLUDE_DIRS})

    #add_package_dependency(ExodusII DEPENDS_ON NetCDF)

   
endif(ExodusII_LIBRARIES AND ExodusII_INCLUDE_DIRS )    

# Define the version
if ( ExodusII_INCLUDE_DIR )
    set(exodus_h ${ExodusII_INCLUDE_DIR}/exodusII.h)
    file(STRINGS ${exodus_h} exodus_version_string REGEX "^#define EX_API_VERS")
    string(REGEX REPLACE "^#define EX_API_VERS ([0-9]+\\.[0-9]+).*$" "\\1" exodus_version "${exodus_version_string}")

    #PRINT_VARIABLE(exodus_version_string)
    #PRINT_VARIABLE(exodus_version)

    set(ExodusII_VERSION "${exodus_version}")

endif()    

# Send useful message if everything is found
find_package_handle_standard_args(ExodusII DEFAULT_MSG
                                           ExodusII_INCLUDE_DIR
                                           ExodusII_LIBRARIES)

# find_package)handle)standard_args should set ExodusII_FOUND but it does not!
if ( ExodusII_LIBRARIES AND ExodusII_INCLUDE_DIRS)
    set(ExodusII_FOUND TRUE)
else()
    set(ExodusII_FOUND FALSE)
endif()


mark_as_advanced(
  ExodusII_VERSION
  ExodusII_INCLUDE_DIR
  ExodusII_INCLUDE_DIRS
  ExodusII_LIBRARY
  ExodusII_LIBRARIES
  ExodusII_LIBRARY_DIR
)
