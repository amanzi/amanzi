# -*- mode: cmake -*-

#
# Amanzi HYPRE Find Module
#
# Usage:
#    Control the search through HYPRE_DIR or setting environment variable
#    HYPRE_ROOT to the HYPRE installation prefix.
#
#    This module does not search default paths! 
#
#    Following variables are set:
#    HYPRE_FOUND            (BOOL)       Flag indicating if HYPRE was found
#    HYPRE_INCLUDE_DIR      (PATH)       Path to the HYPRE include file
#    HYPRE_INCLUDE_DIRS     (LIST)       List of all required include files
#    HYPRE_LIBRARY_DIR      (PATH)       Path to the HYPRE library
#    HYPRE_LIBRARY          (FILE)       HYPRE library
#    HYPRE_LIBRARIES        (LIST)       List of all required HYPRE libraries
#
# #############################################################################

# Standard CMake modules see CMAKE_ROOT/Modules
#include(FindPackageHandleStandardArgs)

# Amanzi CMake functions see <root>/tools/cmake for source
include(PrintVariable)
include(AddPackageDependency)

if (HYPRE_LIBRARIES AND HYPRE_INCLUDE_DIRS)

  # Do nothing. Variables are set. No need to search again

else(HYPRE_LIBRARIES AND HYPRE_INCLUDE_DIRS)

  # Cache variables
  if (HYPRE_DIR)
    set(HYPRE_DIR "${HYPRE_DIR}" CACHE PATH "Path to search for HYPRE include and library files")
  endif()

  if (HYPRE_INCLUDE_DIR)
    set(HYPRE_INCLUDE_DIR "${HYPRE_INCLUDE_DIR}" CACHE PATH "Path to search for HYPRE include files")
  endif()

  if (HYPRE_LIBRARY_DIR)
    set(HYPRE_LIBRARY_DIR "${HYPRE_LIBRARY_DIR}" CACHE PATH "Path to search for HYPRE library files")
  endif()

  set(hypre_inc_names "HYPRE.h")
  if (HYPRE_INCLUDE_DIR)
    # Search for include files in provided HYPRE_INCLUDE_DIR
    if (EXISTS "${HYPRE_INCLUDE_DIR}")

      find_path(hypre_include_path
                NAMES ${hypre_inc_names}
                HINTS ${HYPRE_INCLUDE_DIR}
                NO_DEFAULT_PATH)

      if (NOT hypre_include_path)
        message(SEND_ERROR "Can not locate ${hypre_inc_names} in ${HYPRE_INCLUDE_DIR}")
      endif()

    else()
      message(SEND_ERROR "HYPRE_INCLUDE_DIR=${HYPRE_INCLUDE_DIR} does not exist")
      set(HYPRE_INCLUDE_DIR "HYPRE_INCLUDE_DIR-NOTFOUND")
    endif()

  else() 

    # Search for include files in HYPRE_DIR/include
    set(hypre_inc_suffixes "include")
    if (HYPRE_DIR)
      if (EXISTS "${HYPRE_DIR}")

        find_path(HYPRE_INCLUDE_DIR
                  NAMES ${hypre_inc_names}
                  HINTS ${HYPRE_DIR}
                  PATH_SUFFIXES ${hypre_inc_suffixes}
                  NO_DEFAULT_PATH)

      else()
        message(SEND_ERROR "HYPRE_DIR=${HYPRE_DIR} does not exist")
        set(HYPRE_INCLUDE_DIR "HYPRE_INCLUDE_DIR-NOTFOUND")
      endif()    

    else()

      # Search for include files in default path
      find_path(HYPRE_INCLUDE_DIR
                NAMES ${hypre_inc_names}
                PATH_SUFFIXES ${hypre_inc_suffixes})

    endif()
  endif()

  if (NOT HYPRE_INCLUDE_DIR)
    message(SEND_ERROR "Can not locate HYPRE include directory")
  endif()

  # Search for libraries 
  set(hypre_lib_names "HYPRE")
  if (HYPRE_LIBRARY_DIR)
    # Search for library files in provided HYPRE_LIBRARY_DIR
    if (EXISTS "${HYPRE_LIBRARY_DIR}")

      find_library(_HYPRE_LIBRARY
                   NAMES ${hypre_lib_names}
                   HINTS ${HYPRE_LIBRARY_DIR}
                   NO_DEFAULT_PATH)
    else()
      message(SEND_ERROR "HYPRE_LIBRARY_DIR=${HYPRE_LIBRARY_DIR} does not exist")
      set(HYPRE_LIBRARY "HYPRE_LIBRARY-NOTFOUND")
    endif()

  else() 

    # Search for library files in HYPRE_DIR/lib
    list(APPEND hypre_lib_suffixes "lib" "Lib" "lib64")
    if (HYPRE_DIR)
      if (EXISTS "${HYPRE_DIR}" )

         find_library(_HYPRE_LIBRARY
                      NAMES ${hypre_lib_names}
                      HINTS ${HYPRE_DIR}
                      PATH_SUFFIXES ${hypre_lib_suffixes}
                      NO_DEFAULT_PATH)

      else()
        message(SEND_ERROR "HYPRE_DIR=${HYPRE_DIR} does not exist")
        set(_HYPRE_LIBRARY _HYPRE_LIBRARY-NOTFOUND)
      endif()    

    else()

      # Search for library files in default path
      find_library(_HYPRE_LIBRARY
                   NAMES ${hypre_lib_names}
                   PATH_SUFFIXES ${hypre_lib_suffixes})

    endif()
  endif()

  # Create the target
  if (_HYPRE_LIBRARY)
    set(HYPRE_LIBRARY ${hypre_lib_names})
    get_filename_component(HYPRE_LIBRARY_DIR ${_HYPRE_LIBRARY} DIRECTORY)

    add_imported_library(${HYPRE_LIBRARY}
                         LOCATION ${_HYPRE_LIBRARY}
	                 LINK_LANGUAGES "C")
  else()			   
    message(SEND_ERROR "Can not locate HYPRE library")
  endif()    

  # Update the INCLUDE_DIRS and LIBRARIES variables
  set(HYPRE_INCLUDE_DIRS ${HYPRE_INCLUDE_DIR})
  set(HYPRE_LIBRARIES    ${_HYPRE_LIBRARY})

  # Define the dependent libraries
  set(_HYPRE_DEP_LIBS)

  # -- SuperLU 
  find_library(SuperLU_LIBRARY superlu
               HINTS ${HYPRE_DIR}
               PATH_SUFFIXES lib)

  if (SuperLU_LIBRARY)
    list(APPEND _HYPRE_DEP_LIBS ${SuperLU_LIBRARY})
  else()
    message(SEND_ERROR "Can not locate SuperLU library")
  endif()

  # -- SuperLUDist
  find_library(SuperLUDist_LIBRARY superlu_dist
               HINTS ${HYPRE_DIR}
               PATH_SUFFIXES lib)

  if (SuperLUDist_LIBRARY)
    list(APPEND _HYPRE_DEP_LIBS ${SuperLUDist_LIBRARY})
  else()
    message(SEND_ERROR "Can not locate SuperLUDist library")
  endif()

  # -- ParMetis
  find_library(ParMetis_LIBRARY parmetis
               HINTS ${HYPRE_DIR}
               PATH_SUFFIXES include lib)

  if (ParMetis_LIBRARY)
    list(APPEND _HYPRE_DEP_LIBS ${ParMetis_LIBRARY})
  else()
    message(SEND_ERROR "Can not locate ParMetis library")
  endif()

  # -- METIS
  find_package(METIS QUIET REQUIRED)
  list(APPEND _HYPRE_DEP_LIBS ${METIS_LIBRARY})

  list(APPEND HYPRE_LIBRARIES ${_HYPRE_DEP_LIBS})

  # set target properties   
  set_target_properties(${HYPRE_LIBRARY} PROPERTIES
                       IMPORTED_LINK_INTERFACE_LIBRARIES "${_HYPRE_DEP_LIBS}")

endif(HYPRE_LIBRARIES AND HYPRE_INCLUDE_DIRS)    

# Send useful message if everything is found
 find_package_handle_standard_args(HYPRE DEFAULT_MSG
                                   HYPRE_LIBRARIES
                                   HYPRE_INCLUDE_DIRS)

set(HYPRE_FOUND TRUE)

