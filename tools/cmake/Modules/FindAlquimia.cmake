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

  # Cache variables
  if (Alquimia_DIR)
    set(Alquimia_DIR "${Alquimia_DIR}" CACHE PATH "Path to search for Alquimia include and library files")
  endif()
  set(Alquimia_INCLUDE_DIR ${Alquimia_DIR}/include)
  set(Alquimia_LIBRARY_DIR ${Alquimia_DIR}/lib)
  set(Alquimia_TARGET alquimia)

  if (Alquimia_INCLUDE_DIR)
    set(Alquimia_INCLUDE_DIR "${Alquimia_INCLUDE_DIR}" CACHE PATH "Path to search for Alquimia include files")
  endif()
  find_library(_Alquimia_LIBRARY
               NAMES alquimia
               PATHS ${Alquimia_LIBRARY_DIR})
  if (Alquimia_LIBRARY_DIR)
    set(Alquimia_LIBRARY_DIR "${Alquimia_LIBRARY_DIR}" CACHE PATH "Path to search for Alquimia library files")
  endif()
  if (_Alquimia_LIBRARY)
    add_imported_library(${Alquimia_TARGET} 
                         LOCATION ${_Alquimia_LIBRARY}
                         LINK_LANGUAGES "C")
    set(Alquimia_LIBRARY ${Alquimia_TARGET})
  endif()    
    
  # Search for include files
  # Search order preference:
  #  (1) Alquimia_INCLUDE_DIR - check existence of path AND if the include files exist
  #  (2) Alquimia_DIR/<include>
  #  (3) Default CMake paths See cmake --html-help=out.html file for more information.
  #
  set(alquimia_inc_names "alquimia.h")
  if (Alquimia_INCLUDE_DIR)
    if (EXISTS "${Alquimia_INCLUDE_DIR}")

      find_path(alquimia_test_include_path
                NAMES ${alquimia_inc_names}
                HINTS ${Alquimia_INCLUDE_DIR}
                NO_DEFAULT_PATH)

      if (NOT alquimia_test_include_path)
        message(SEND_ERROR "Can not locate ${alquimia_inc_names} in ${Alquimia_INCLUDE_DIR}")
      endif()
      set(Alquimia_INCLUDE_DIR "${alquimia_test_include_path}")

    else()
      message(SEND_ERROR "Alquimia_INCLUDE_DIR=${Alquimia_INCLUDE_DIR} does not exist")
      set(Alquimia_INCLUDE_DIR "Alquimia_INCLUDE_DIR-NOTFOUND")
    endif()

  else() 

    set(alquimia_inc_suffixes "include/alquimia")
    if (Alquimia_DIR)
      if (EXISTS "${Alquimia_DIR}")

        find_path(Alquimia_INCLUDE_DIR
                  NAMES ${alquimia_inc_names}
                  HINTS ${Alquimia_DIR}
                  PATH_SUFFIXES ${alquimia_inc_suffixes}
                  NO_DEFAULT_PATH)

      else()
        message(SEND_ERROR "Alquimia_DIR=${Alquimia_DIR} does not exist")
        set(Alquimia_INCLUDE_DIR "Alquimia_INCLUDE_DIR-NOTFOUND")
      endif()    

    else()

      find_path(Alquimia_INCLUDE_DIR
                NAMES ${alquimia_inc_names}
                PATH_SUFFIXES ${alquimia_inc_suffixes})

    endif()
  endif()

  if (NOT Alquimia_INCLUDE_DIR)
    message(SEND_ERROR "Can not locate Alquimia include directory")
  endif()

  # Search for libraries 
  # Search order preference:
  #  (1) Alquimia_LIBRARY_DIR - check existence of path AND if the library file exists
  #  (2) Alquimia_DIR/<lib,Lib>
  #  (3) Default CMake paths See cmake --html-help=out.html file for more information.
  #
  set(alquimia_lib_names "alquimia")
  if (Alquimia_LIBRARY_DIR)
    if (EXISTS "${Alquimia_LIBRARY_DIR}")

      find_library(_Alquimia_LIBRARY
                   NAMES ${alquimia_lib_names}
                   HINTS ${Alquimia_LIBRARY_DIR}
                   NO_DEFAULT_PATH)

    else()
      message(SEND_ERROR "Alquimia_LIBRARY_DIR=${Alquimia_LIBRARY_DIR} does not exist")
      set(_Alquimia_LIBRARY "Alquimia_LIBRARY-NOTFOUND")
      set(_Alquimia_Fortran_LIBRARY "Alquimia_Fortran_LIBRARY-NOTFOUND")
    endif()

  else() 

    list(APPEND alquimia_lib_suffixes "lib" "Lib")
    if (Alquimia_DIR)
      if (EXISTS "${Alquimia_DIR}")

        find_library(_Alquimia_LIBRARY
                     NAMES ${alquimia_lib_names}
                     HINTS ${Alquimia_DIR}
                     PATH_SUFFIXES ${alquimia_lib_suffixes}
                     NO_DEFAULT_PATH)
                
      else()
        message(SEND_ERROR "Alquimia_DIR=${Alquimia_DIR} does not exist")
        set(Alquimia_LIBRARY "Alquimia_LIBRARY-NOTFOUND")
        set(Alquimia_Fortran_LIBRARY "Alquimia_Fortran_LIBRARY-NOTFOUND")
      endif()    

    else()

      find_library(_Alquimia_LIBRARY
                   NAMES ${alquimia_lib_names}
                   PATH_SUFFIXES ${alquimia_lib_suffixes})

    endif()
  endif()

  # Create the library target store the name in Alquimia_LIBRARY
  if ( _Alquimia_LIBRARY )
    set(Alquimia_LIBRARY alquimia)
    add_imported_library(${Alquimia_LIBRARY}
                        LOCATION ${Alquimia_LIBRARY})
  else()
    message(SEND_ERROR "Can not locate Alquimia library")
  endif()

  # Define prerequisite packages
  set(Alquimia_INCLUDE_DIRS ${Alquimia_INCLUDE_DIR})
  set(Alquimia_LIBRARIES    ${Alquimia_LIBRARY})

  # Alquimia generates a CMake configuration file that contains the
  # required TPLs. I use an include here instead of find_package
  # to prevent a recursive call.
  if (Alquimia_DIR)
    set(Alquimia_CMAKE_CONFIG_FILE ${Alquimia_DIR}/share/alquimia/alquimia.cmake)
    if (EXISTS ${Alquimia_CMAKE_CONFIG_FILE})
      include(${Alquimia_CMAKE_CONFIG_FILE})
  # Define the LIBRARIES and INCLUDE_DIRS
  set(Alquimia_INCLUDE_DIRS ${Alquimia_INCLUDE_DIR})
  set(Alquimia_LIBRARIES ${Alquimia_LIBRARY} ${CrunchTope_LIBRARY} ${PFLOTRAN_LIBRARIES})

      # Include paths
      if (Alquimia_PACKAGE_INCLUDES)
	      list(APPEND Alquimia_INCLUDE_DIRS ${Alquimia_PACKAGE_INCLUDES})
	      list(REMOVE_DUPLICATES Alquimia_INCLUDE_DIRS)
      endif()

      # TPL libraries, some of the items in this list are not defined!
      if (Alquimia_PACKAGE_LIBS)
	      foreach(lib ${Alquimia_PACKAGE_LIBS})
	        if (lib)
	          list(APPEND Alquimia_LIBRARIES ${lib})
	        endif()
        endforeach()
      endif()  
    endif()
  endif()  
   
endif(Alquimia_LIBRARIES AND Alquimia_INCLUDE_DIRS) 

# Send useful message if everything is found
find_package_handle_standard_args(Alquimia DEFAULT_MSG
                                  Alquimia_INCLUDE_DIR
                                  Alquimia_LIBRARIES)

# find_package)handle)standard_args should set Alquimia_FOUND but it does not!
if (Alquimia_LIBRARIES AND Alquimia_INCLUDE_DIRS)
  set(Alquimia_FOUND TRUE)
else()
  set(Alquimia_FOUND FALSE)
endif()

mark_as_advanced(
  Alquimia_INCLUDE_DIR
  Alquimia_INCLUDE_DIRS
  Alquimia_LIBRARY
  Alquimia_LIBRARIES
  Alquimia_LIBRARY_DIR
)
