# -*- mode: cmake -*-

#
# Amanzi Alquimia Find Module
#
# Usage:
#    Control the search through Alquimia_DIR or setting environment variable
#    Alquimia_ROOT to the Alquimia installation prefix.
#
#    This module does not search default paths! 
#
#    Following variables are set:
#    Alquimia_FOUND            (BOOL)       Flag indicating if Alquimia was found
#    Alquimia_INCLUDE_DIR      (PATH)       Path to the Alquimia include file
#    Alquimia_INCLUDE_DIRS     (LIST)       List of all required include files
#    Alquimia_LIBRARY_DIR      (PATH)       Path to the Alquimia library
#    Alquimia_LIBRARIES        (LIST)       List of all required Alquimia libraries
#
# #############################################################################

include(AddImportedLibrary)

# Amanzi CMake functions see <root>/tools/cmake for source
include(PrintVariable)
include(AddPackageDependency)

if (Alquimia_LIBRARIES AND Alquimia_INCLUDE_DIRS)

  # Do nothing. Variables are set. No need to search again

elseif(Alquimia_DIR)

  set(Alquimia_INCLUDE_DIR ${Alquimia_DIR}/include)
  set(Alquimia_LIBRARY_DIR ${Alquimia_DIR}/lib)
  set(Alquimia_TARGET alquimia)
  find_library(_Alquimia_LIBRARY
               NAMES alquimia
               PATHS ${Alquimia_LIBRARY_DIR})

  if (_Alquimia_LIBRARY)
    add_imported_library(${Alquimia_TARGET} 
                         LOCATION ${_Alquimia_LIBRARY}
                         LINK_LANGUAGES "C")
    set(Alquimia_LIBRARY ${Alquimia_TARGET})
  endif()    
  
  # Define the LIBRARIES and INCLUDE_DIRS
  set(Alquimia_INCLUDE_DIRS ${Alquimia_INCLUDE_DIR})
  set(Alquimia_LIBRARIES ${Alquimia_LIBRARY} ${CrunchTope_LIBRARY} ${PFLOTRAN_LIBRARIES})
  
  # Search for include files
  # Search order preference:
  #  (1) ALQUIMIA_INCLUDE_DIR - check existence of path AND if the include files exist
  #  (2) ALQUIMIA_DIR/<include>
  #  (3) Default CMake paths See cmake --html-help=out.html file for more information.
  #
  set(alquimia_inc_names "alquimia.h")
  if (ALQUIMIA_INCLUDE_DIR)
    if (EXISTS "${ALQUIMIA_INCLUDE_DIR}")

      find_path(alquimia_test_include_path
                NAMES ${alquimia_inc_names}
                HINTS ${ALQUIMIA_INCLUDE_DIR}
                NO_DEFAULT_PATH)

      if (NOT alquimia_test_include_path)
        message(SEND_ERROR "Can not locate ${alquimia_inc_names} in ${ALQUIMIA_INCLUDE_DIR}")
      endif()
      set(ALQUIMIA_INCLUDE_DIR "${alquimia_test_include_path}")

    else()
      message(SEND_ERROR "ALQUIMIA_INCLUDE_DIR=${ALQUIMIA_INCLUDE_DIR} does not exist")
      set(ALQUIMIA_INCLUDE_DIR "ALQUIMIA_INCLUDE_DIR-NOTFOUND")
    endif()

  else() 

    set(alquimia_inc_suffixes "include/alquimia")
    if (ALQUIMIA_DIR)
      if (EXISTS "${ALQUIMIA_DIR}")

        find_path(ALQUIMIA_INCLUDE_DIR
                  NAMES ${alquimia_inc_names}
                  HINTS ${ALQUIMIA_DIR}
                  PATH_SUFFIXES ${alquimia_inc_suffixes}
                  NO_DEFAULT_PATH)

      else()
        message(SEND_ERROR "ALQUIMIA_DIR=${ALQUIMIA_DIR} does not exist")
        set(ALQUIMIA_INCLUDE_DIR "ALQUIMIA_INCLUDE_DIR-NOTFOUND")
      endif()    

    else()

      find_path(ALQUIMIA_INCLUDE_DIR
                NAMES ${alquimia_inc_names}
                PATH_SUFFIXES ${alquimia_inc_suffixes})

    endif()
  endif()

  if (NOT ALQUIMIA_INCLUDE_DIR)
    message(SEND_ERROR "Can not locate Alquimia include directory")
  endif()

  # Search for libraries 
  # Search order preference:
  #  (1) ALQUIMIA_LIBRARY_DIR - check existence of path AND if the library file exists
  #  (2) ALQUIMIA_DIR/<lib,Lib>
  #  (3) Default CMake paths See cmake --html-help=out.html file for more information.
  #
  set(alquimia_lib_names "alquimia")
  if (ALQUIMIA_LIBRARY_DIR)
    if (EXISTS "${ALQUIMIA_LIBRARY_DIR}")

      find_library(_ALQUIMIA_LIBRARY
                   NAMES ${alquimia_lib_names}
                   HINTS ${ALQUIMIA_LIBRARY_DIR}
                   NO_DEFAULT_PATH)

    else()
      message(SEND_ERROR "ALQUIMIA_LIBRARY_DIR=${ALQUIMIA_LIBRARY_DIR} does not exist")
      set(_ALQUIMIA_LIBRARY "ALQUIMIA_LIBRARY-NOTFOUND")
      set(_ALQUIMIA_Fortran_LIBRARY "ALQUIMIA_Fortran_LIBRARY-NOTFOUND")
    endif()

  else() 

    list(APPEND alquimia_lib_suffixes "lib" "Lib")
    if (ALQUIMIA_DIR)
      if (EXISTS "${ALQUIMIA_DIR}")

        find_library(_ALQUIMIA_LIBRARY
                     NAMES ${alquimia_lib_names}
                     HINTS ${ALQUIMIA_DIR}
                     PATH_SUFFIXES ${alquimia_lib_suffixes}
                     NO_DEFAULT_PATH)
                
      else()
        message(SEND_ERROR "ALQUIMIA_DIR=${ALQUIMIA_DIR} does not exist")
        set(ALQUIMIA_LIBRARY "ALQUIMIA_LIBRARY-NOTFOUND")
        set(ALQUIMIA_Fortran_LIBRARY "ALQUIMIA_Fortran_LIBRARY-NOTFOUND")
      endif()    

    else()

      find_library(_ALQUIMIA_LIBRARY
                   NAMES ${alquimia_lib_names}
                   PATH_SUFFIXES ${alquimia_lib_suffixes})

    endif()
  endif()

  # Create the library target store the name in ALQUIMIA_LIBRARY
  if ( _ALQUIMIA_LIBRARY )
    set(ALQUIMIA_LIBRARY alquimia)
    add_imported_library(${ALQUIMIA_LIBRARY}
                        LOCATION ${ALQUIMIA_LIBRARY})
  else()
    message(SEND_ERROR "Can not locate ALQUIMIA library")
  endif()

  # Define prerequisite packages
  set(ALQUIMIA_INCLUDE_DIRS ${ALQUIMIA_INCLUDE_DIR})
  set(ALQUIMIA_LIBRARIES    ${ALQUIMIA_LIBRARY})

  # Alquimia generates a CMake configuration file that contains the
  # required TPLs. I use an include here instead of find_package
  # to prevent a recursive call.
  if (ALQUIMIA_DIR)
    set(ALQUIMIA_CMAKE_CONFIG_FILE ${ALQUIMIA_DIR}/share/alquimia/alquimia.cmake)
    if (EXISTS ${ALQUIMIA_CMAKE_CONFIG_FILE})
      include(${ALQUIMIA_CMAKE_CONFIG_FILE})
endif(Alquimia_LIBRARIES AND Alquimia_INCLUDE_DIRS)    
   

# Send useful message if everything is found
find_package_handle_standard_args(Alquimia DEFAULT_MSG
                                  Alquimia_INCLUDE_DIRS
                                  Alquimia_LIBRARIES)

# find_package)handle)standard_args should set Alquimia_FOUND but it does not!
if (Alquimia_LIBRARIES AND Alquimia_INCLUDE_DIRS)
  set(Alquimia_FOUND TRUE)
else()
  set(Alquimia_FOUND FALSE)
endif()

mark_as_advanced(
  Alquimia_INCLUDE_DIR
  Alquimia_LIBRARY_DIR
  Alquimia_LIBRARY
  Alquimia_LIBRARIES
)
