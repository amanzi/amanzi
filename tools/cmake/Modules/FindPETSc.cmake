# -*- mode: cmake -*-

#
# Amanzi PETSc Find Module
#
# Usage:
#    Control the search through PETSc_DIR or setting environment variable
#    PETSc_ROOT to the PETSc installation prefix.
#
#    This module does not search default paths! 
#
#    Following variables are set:
#    PETSc_FOUND            (BOOL)       Flag indicating if PETSc was found
#    PETSc_INCLUDE_DIR      (PATH)       Path to the PETSc include file
#    PETSc_INCLUDE_DIRS     (LIST)       List of all required include files
#    PETSc_LIBRARY_DIR      (PATH)       Path to the PETSc library
#    PETSc_LIBRARY          (FILE)       PETSc library
#    PETSc_LIBRARIES        (LIST)       List of all required PETSc libraries
#
#    Additional variables
#    PETSc_VERSION          (STRING)     PETSc Version string
#
# #############################################################################

# Standard CMake modules see CMAKE_ROOT/Modules
include(FindPackageHandleStandardArgs)

include(AddImportedLibrary)

# Amanzi CMake functions see <root>/tools/cmake for source
include(PrintVariable)
include(AddPackageDependency)

if (PETSc_LIBRARIES AND PETSc_INCLUDE_DIRS)

  # Do nothing. Variables are set. No need to search again

else(PETSc_LIBRARIES AND PETSc_INCLUDE_DIRS)

  # Cache variables
  if (PETSc_DIR)
    set(PETSc_DIR "${PETSc_DIR}" CACHE PATH "Path to search for PETSc include and library files")
  endif()

  if (PETSc_INCLUDE_DIR)
    set(PETSc_INCLUDE_DIR "${PETSc_INCLUDE_DIR}" CACHE PATH "Path to search for PETSc include files")
  endif()

  if (PETSc_LIBRARY_DIR)
    set(PETSc_LIBRARY_DIR "${PETSc_LIBRARY_DIR}" CACHE PATH "Path to search for PETSc library files")
  endif()

    
  # Search for include files
  # Search order preference:
  #  (1) PETSc_INCLUDE_DIR - check existence of path AND if the include files exist
  #  (2) PETSc_DIR/<include>
  #  (3) Default CMake paths See cmake --html-help=out.html file for more information.
  #
  set(petsc_inc_names "petsc.h")
  if (PETSc_INCLUDE_DIR)
    if (EXISTS "${PETSc_INCLUDE_DIR}")

      find_path(petsc_test_include_path
                NAMES ${petsc_inc_names}
                HINTS ${PETSc_INCLUDE_DIR}
                NO_DEFAULT_PATH)

      if (NOT petsc_test_include_path)
        message(SEND_ERROR "Can not locate ${petsc_inc_names} in ${PETSc_INCLUDE_DIR}")
      endif()
      set(PETSc_INCLUDE_DIR "${petsc_test_include_path}")

    else()
      message(SEND_ERROR "PETSc_INCLUDE_DIR=${PETSc_INCLUDE_DIR} does not exist")
      set(PETSc_INCLUDE_DIR "PETSc_INCLUDE_DIR-NOTFOUND")
    endif()

  else() 

    set(petsc_inc_suffixes "include")
    if (PETSc_DIR)
      if (EXISTS "${PETSc_DIR}")

        find_path(PETSc_INCLUDE_DIR
                  NAMES ${petsc_inc_names}
                  HINTS ${PETSc_DIR}
                  PATH_SUFFIXES ${petsc_inc_suffixes}
                  NO_DEFAULT_PATH)
                
      else()
        message(SEND_ERROR "PETSc_DIR=${PETSc_DIR} does not exist")
        set(PETSc_INCLUDE_DIR "PETSc_INCLUDE_DIR-NOTFOUND")
      endif()    

    else()

      find_path(PETSc_INCLUDE_DIR
                NAMES ${petsc_inc_names}
                PATH_SUFFIXES ${petsc_inc_suffixes})

    endif()
  endif()

  if (NOT PETSc_INCLUDE_DIR)
    message(SEND_ERROR "Can not locate PETSc include directory")
  endif()

  # Search for libraries 
  # Search order preference:
  #  (1) PETSc_LIBRARY_DIR - check existence of path AND if the library file exists
  #  (2) PETSc_DIR/<lib,Lib>
  #  (3) Default CMake paths See cmake --html-help=out.html file for more information.
  #
  set(petsc_lib_names "petsc")
  if (PETSc_LIBRARY_DIR)
    if (EXISTS "${PETSc_LIBRARY_DIR}")

      find_library(_PETSc_LIBRARY
                   NAMES ${petsc_lib_names}
                   HINTS ${PETSc_LIBRARY_DIR}
                   NO_DEFAULT_PATH)

    else()
      message(SEND_ERROR "PETSc_LIBRARY_DIR=${PETSc_LIBRARY_DIR} does not exist")
      set(_PETSc_LIBRARY "PETSc_LIBRARY-NOTFOUND")
      set(_PETSc_Fortran_LIBRARY "PETSc_Fortran_LIBRARY-NOTFOUND")
    endif()

  else() 

    list(APPEND petsc_lib_suffixes "lib" "Lib")
    if (PETSc_DIR)
      if (EXISTS "${PETSc_DIR}")

        find_library(_PETSc_LIBRARY
                     NAMES ${petsc_lib_names}
                     HINTS ${PETSc_DIR}
                     PATH_SUFFIXES ${petsc_lib_suffixes}
                     NO_DEFAULT_PATH)
                
      else()
        message(SEND_ERROR "PETSc_DIR=${PETSc_DIR} does not exist")
        set(PETSc_LIBRARY "PETSc_LIBRARY-NOTFOUND")
        set(PETSc_Fortran_LIBRARY "PETSc_Fortran_LIBRARY-NOTFOUND")
      endif()    

    else()

      find_library(_PETSc_LIBRARY
                   NAMES ${petsc_lib_names}
                   PATH_SUFFIXES ${petsc_lib_suffixes})

    endif()
  endif()

  # Create the library target store the name in PETSc_LIBRARY
  if ( _PETSc_LIBRARY )
    set(PETSc_LIBRARY petsc)
    add_imported_library(${PETSc_LIBRARY}
                        LOCATION ${_PETSc_LIBRARY})
  else()
    message(SEND_ERROR "Can not locate PETSc library")
  endif()

  # Define prerequisite packages
  set(PETSc_INCLUDE_DIRS ${PETSc_INCLUDE_DIR})
  set(PETSc_LIBRARIES    ${PETSc_LIBRARY})

  # PETSc generates a CMake configuration file that contains the
  # required TPLs. I use an include here instead of find_package
  # to prevent a recursive call.
  if (PETSc_DIR)
    set(PETSc_CMAKE_CONFIG_FILE ${PETSc_DIR}/lib/petsc/conf/PETScBuildInternal.cmake)
    if (EXISTS ${PETSc_CMAKE_CONFIG_FILE})
      include(${PETSc_CMAKE_CONFIG_FILE})

      # Include paths
      if (PETSC_PACKAGE_INCLUDES)
        list(APPEND PETSc_INCLUDE_DIRS ${PETSC_PACKAGE_INCLUDES})
        list(REMOVE_DUPLICATES PETSc_INCLUDE_DIRS)
      endif()

      # TPL libraries, some of the items in this list are not defined!
      if (PETSC_PACKAGE_LIBS)
	    foreach(lib ${PETSC_PACKAGE_LIBS})
	      if (lib)
	        list(APPEND PETSc_LIBRARIES ${lib})
	      endif()
        endforeach()
      endif()  

      # set target properties   
      set(petsc_deps HYPRE SUPERLU_DIST SUPERLU PARMETIS METIS)
      foreach(lib ${petsc_deps})
        message(STATUS ">>>>> JDM: ${lib} ${PETSC_HAVE_${lib}}")
        if (PETSC_HAVE_${lib})
          set_target_properties(${PETSc_LIBRARY} PROPERTIES
                                INTERFACE_LINK_LIBRARIES "${PETSC_${lib}_LIB}")
        endif()
      endforeach()
    endif()
  endif()  
   
endif(PETSc_LIBRARIES AND PETSc_INCLUDE_DIRS)

# Define the version
if (NOT PETSc_VERSION)
  set(PETSc_VERSION "")
  if (PETSc_INCLUDE_DIR)
    set(petscversion_h ${PETSc_INCLUDE_DIR}/petscversion.h)
    if (EXISTS ${petscversion_h})
      set(version_labels MAJOR MINOR PATCH)
      foreach(label ${version_labels})
        set(regexp_target "\#define PETSC_VERSION_${label}[ \t]+")
        file(STRINGS ${petscversion_h} version_string REGEX "^${regexp_target}")
        string(REGEX REPLACE "${regexp_target}\([0-9]+\)[ \t\n\r]*" "\\1" ver_num ${version_string})

        if (ver_num)
	  if (PETSc_VERSION)
	    set(PETSc_VERSION "${PETSc_VERSION}.${ver_num}")
          else()
    	    set(PETSc_VERSION ${ver_num})
	  endif()
        endif()  
      endforeach()  
    endif()
  endif()

endif()    

# Send useful message if everything is found
find_package_handle_standard_args(PETSc DEFAULT_MSG
                                  PETSc_INCLUDE_DIR
                                  PETSc_LIBRARIES)

mark_as_advanced(
  PETSc_INCLUDE_DIR
  PETSc_INCLUDE_DIRS
  PETSc_LIBRARY
  PETSc_LIBRARIES
  PETSc_LIBRARY_DIR
)
