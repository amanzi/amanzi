# -*- mode: cmake -*-

#
# Amanzi PETSC Find Module
#
# Usage:
#    Control the search through PETSC_DIR or setting environment variable
#    PETSC_ROOT to the PETSC installation prefix.
#
#    This module does not search default paths! 
#
#    Following variables are set:
#    PETSC_FOUND            (BOOL)       Flag indicating if PETSC was found
#    PETSC_INCLUDE_DIR      (PATH)       Path to the PETSC include file
#    PETSC_INCLUDE_DIRS     (LIST)       List of all required include files
#    PETSC_LIBRARY_DIR      (PATH)       Path to the PETSC library
#    PETSC_LIBRARY          (FILE)       PETSC library
#    PETSC_LIBRARIES        (LIST)       List of all required PETSC libraries
#
#    Additional variables
#    PETSC_VERSION          (STRING)     PETSC Version string
#
# #############################################################################

# Standard CMake modules see CMAKE_ROOT/Modules
include(FindPackageHandleStandardArgs)

include(AddImportedLibrary)

# Amanzi CMake functions see <root>/tools/cmake for source
include(PrintVariable)
include(AddPackageDependency)

if ( PETSC_LIBRARIES AND PETSC_INCLUDE_DIRS )

    # Do nothing. Variables are set. No need to search again

else(PETSC_LIBRARIES AND PETSC_INCLUDE_DIRS)

    # Cache variables
    if(PETSC_DIR)
        set(PETSC_DIR "${PETSC_DIR}" CACHE PATH "Path to search for PETSC include and library files")
    endif()

    if(PETSC_INCLUDE_DIR)
        set(PETSC_INCLUDE_DIR "${PETSC_INCLUDE_DIR}" CACHE PATH "Path to search for PETSC include files")
    endif()

    if(PETSC_LIBRARY_DIR)
        set(PETSC_LIBRARY_DIR "${PETSC_LIBRARY_DIR}" CACHE PATH "Path to search for PETSC library files")
    endif()

    
    # Search for include files
    # Search order preference:
    #  (1) PETSC_INCLUDE_DIR - check existence of path AND if the include files exist
    #  (2) PETSC_DIR/<include>
    #  (3) Default CMake paths See cmake --html-help=out.html file for more information.
    #
    set(petsc_inc_names "petsc.h")
    if (PETSC_INCLUDE_DIR)

        if (EXISTS "${PETSC_INCLUDE_DIR}")

            find_path(petsc_test_include_path
                      NAMES ${petsc_inc_names}
                      HINTS ${PETSC_INCLUDE_DIR}
                      NO_DEFAULT_PATH)
            if(NOT petsc_test_include_path)
                message(SEND_ERROR "Can not locate ${petsc_inc_names} in ${PETSC_INCLUDE_DIR}")
            endif()
            set(PETSC_INCLUDE_DIR "${petsc_test_include_path}")

        else()
            message(SEND_ERROR "PETSC_INCLUDE_DIR=${PETSC_INCLUDE_DIR} does not exist")
            set(PETSC_INCLUDE_DIR "PETSC_INCLUDE_DIR-NOTFOUND")
        endif()

   else() 

        set(petsc_inc_suffixes "include")
        if(PETSC_DIR)

            if (EXISTS "${PETSC_DIR}" )

                find_path(PETSC_INCLUDE_DIR
                          NAMES ${petsc_inc_names}
                          HINTS ${PETSC_DIR}
                          PATH_SUFFIXES ${petsc_inc_suffixes}
                          NO_DEFAULT_PATH)

            else()
                 message(SEND_ERROR "PETSC_DIR=${PETSC_DIR} does not exist")
                 set(PETSC_INCLUDE_DIR "PETSC_INCLUDE_DIR-NOTFOUND")
            endif()    


        else()

            find_path(PETSC_INCLUDE_DIR
                      NAMES ${petsc_inc_names}
                      PATH_SUFFIXES ${petsc_inc_suffixes})

        endif()

    endif()

    if ( NOT PETSC_INCLUDE_DIR )
        message(SEND_ERROR "Can not locate PETSC include directory")
    endif()

    # Search for libraries 
    # Search order preference:
    #  (1) PETSC_LIBRARY_DIR - check existence of path AND if the library file exists
    #  (2) PETSC_DIR/<lib,Lib>
    #  (3) Default CMake paths See cmake --html-help=out.html file for more information.
    #
    set(petsc_lib_names "petsc")
    if (PETSC_LIBRARY_DIR)

        if (EXISTS "${PETSC_LIBRARY_DIR}")

            find_library(_PETSC_LIBRARY
                         NAMES ${petsc_lib_names}
                         HINTS ${PETSC_LIBRARY_DIR}
                         NO_DEFAULT_PATH)

        else()
            message(SEND_ERROR "PETSC_LIBRARY_DIR=${PETSC_LIBRARY_DIR} does not exist")
            set(_PETSC_LIBRARY "PETSC_LIBRARY-NOTFOUND")
            set(_PETSC_Fortran_LIBRARY "PETSC_Fortran_LIBRARY-NOTFOUND")
        endif()

    else() 

        list(APPEND petsc_lib_suffixes "lib" "Lib")
        if(PETSC_DIR)

            if (EXISTS "${PETSC_DIR}" )

                find_library(_PETSC_LIBRARY
                             NAMES ${petsc_lib_names}
                             HINTS ${PETSC_DIR}
                             PATH_SUFFIXES ${petsc_lib_suffixes}
                             NO_DEFAULT_PATH)
                
            else()
                 message(SEND_ERROR "PETSC_DIR=${PETSC_DIR} does not exist")
                 set(PETSC_LIBRARY "PETSC_LIBRARY-NOTFOUND")
                 set(PETSC_Fortran_LIBRARY "PETSC_Fortran_LIBRARY-NOTFOUND")
            endif()    


        else()

            find_library(_PETSC_LIBRARY
                         NAMES ${petsc_lib_names}
                         PATH_SUFFIXES ${petsc_lib_suffixes})

        endif()

    endif()

    # Create the library target store the name in PETSC_LIBRARY
    if ( _PETSC_LIBRARY )
        set(PETSC_LIBRARY petsc)
        add_imported_library(${PETSC_LIBRARY}
                     LOCATION ${_PETSC_LIBRARY})
    else()
        message(SEND_ERROR "Can not locate PETSC library")
    endif()

    # Define prerequisite packages
    set(PETSC_INCLUDE_DIRS ${PETSC_INCLUDE_DIR})
    set(PETSC_LIBRARIES    ${PETSC_LIBRARY})

    # PETSc generates a CMake configuration file that contains the
    # required TPLs. I use an include here instead of find_package
    # to prevent a recursive call.
    if (PETSC_DIR)
      set(PETSC_CMAKE_CONFIG_FILE ${PETSC_DIR}/conf/PETScConfig.cmake)
      if ( EXISTS ${PETSC_CMAKE_CONFIG_FILE})
	include(${PETSC_CMAKE_CONFIG_FILE})

	# Include paths
	if(PETSC_PACKAGE_INCLUDES)
	  list(APPEND PETSC_INCLUDE_DIRS ${PETSC_PACKAGE_INCLUDES})
	  list(REMOVE_DUPLICATES PETSC_INCLUDE_DIRS)
	endif()

	# TPL libraries, some of the items in this list are not defined!
	if ( PETSC_PACKAGE_LIBS)
	  foreach(lib ${PETSC_PACKAGE_LIBS})
	    if(lib)
	      list(APPEND PETSC_LIBRARIES ${lib})
	    endif()
	  endforeach()
	endif()  

      endif()

    endif()  
   
endif(PETSC_LIBRARIES AND PETSC_INCLUDE_DIRS )    

# Define the version
if ( NOT PETSC_VERSION )
  set(PETSC_VERSION "")
  if ( PETSC_INCLUDE_DIR )
      set(petscversion_h ${PETSC_INCLUDE_DIR}/petscversion.h)
      if (EXISTS ${petscversion_h})
        set(version_labels MAJOR MINOR PATCH)
        foreach(label ${version_labels})
	  set(regexp_target "\#define PETSC_VERSION_${label}[ \t]+")
	  file(STRINGS ${petscversion_h} version_string REGEX "^${regexp_target}")
	  #print_variable(version_string)
	  string(REGEX REPLACE "${regexp_target}\([0-9]+\)[ \t\n\r]*" "\\1" ver_num ${version_string})
	  #print_variable(ver_num)
	  if(ver_num)
	    if (PETSC_VERSION)
	      set(PETSC_VERSION "${PETSC_VERSION}.${ver_num}")
            else()
    	      set(PETSC_VERSION ${ver_num})
	    endif()
	  endif()  
	endforeach()  
      endif()
 endif()

 #print_variable(PETSC_VERSION)

endif()    

# Send useful message if everything is found
find_package_handle_standard_args(PETSC DEFAULT_MSG
                                        PETSC_INCLUDE_DIR
                                        PETSC_LIBRARIES)

mark_as_advanced(
  PETSC_INCLUDE_DIR
  PETSC_INCLUDE_DIRS
  PETSC_LIBRARY
  PETSC_LIBRARIES
  PETSC_LIBRARY_DIR
)
