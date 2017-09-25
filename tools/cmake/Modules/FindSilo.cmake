# -*- mode: cmake -*-

#
# Amanzi Silo Find Module
#
# Usage:
#    Control the search through Silo_DIR or setting environment variable
#    Silo_ROOT to the Silo installation prefix.
#
#    This module does not search default paths! 
#
#    Following variables are set:
#    Silo_FOUND            (BOOL)       Flag indicating if Silo was found
#    Silo_INCLUDE_DIR      (PATH)       Path to the Silo include file
#    Silo_INCLUDE_DIRS     (LIST)       List of all required include files
#    Silo_LIBRARY_DIR      (PATH)       Path to the Silo library
#    Silo_LIBRARY          (FILE)       Silo library
#    Silo_LIBRARIES        (LIST)       List of all required Silo libraries
#
# #############################################################################

# Standard CMake modules see CMAKE_ROOT/Modules
#include(FindPackageHandleStandardArgs)

# Amanzi CMake functions see <root>/tools/cmake for source
include(PrintVariable)
include(AddPackageDependency)

if (Silo_LIBRARIES AND Silo_INCLUDE_DIRS)

  # Do nothing. Variables are set. No need to search again

else(Silo_LIBRARIES AND Silo_INCLUDE_DIRS)

  # Cache variables
  if (Silo_DIR)
    set(Silo_DIR "${Silo_DIR}" CACHE PATH "Path to search for Silo include and library files")
  endif()

  if (Silo_INCLUDE_DIR)
    set(Silo_INCLUDE_DIR "${Silo_INCLUDE_DIR}" CACHE PATH "Path to search for Silo include files")
  endif()

  if (Silo_LIBRARY_DIR)
    set(Silo_LIBRARY_DIR "${Silo_LIBRARY_DIR}" CACHE PATH "Path to search for Silo library files")
  endif()

  # Search for include files
  # Search order preference:
  #  (1) Silo_INCLUDE_DIR - check existence of path AND if the include files exist
  #  (2) Silo_DIR/<include>
  #  (3) Default CMake paths See cmake --html-help=out.html file for more information.
  #
  set(silo_inc_names "silo.h")
  if (Silo_INCLUDE_DIR)
    if (EXISTS "${Silo_INCLUDE_DIR}")

      find_path(test_silo_include_path
                NAMES ${silo_inc_names}
                HINTS ${Silo_INCLUDE_DIR}
                NO_DEFAULT_PATH)
      if (NOT test_silo_include_path)
         message(SEND_ERROR "Can not locate ${silo_inc_names} in ${Silo_INCLUDE_DIR}")
      endif()
      set(Silo_INCLUDE_DIR "${test_silo_include_path}")

    else()
       message(SEND_ERROR "Silo_INCLUDE_DIR=${Silo_INCLUDE_DIR} does not exist")
       set(Silo_INCLUDE_DIR "Silo_INCLUDE_DIR-NOTFOUND")
    endif()

  else() 

    set(silo_inc_suffixes "include")
    if (Silo_DIR)
      if (EXISTS "${Silo_DIR}")

        find_path(Silo_INCLUDE_DIR
                  NAMES ${silo_inc_names}
                  HINTS ${Silo_DIR}
                  PATH_SUFFIXES ${silo_inc_suffixes}
                  NO_DEFAULT_PATH)

      else()
        message(SEND_ERROR "Silo_DIR=${Silo_DIR} does not exist")
        set(Silo_INCLUDE_DIR "Silo_INCLUDE_DIR-NOTFOUND")
      endif()    

    else()

      find_path(Silo_INCLUDE_DIR
                NAMES ${silo_inc_names}
                PATH_SUFFIXES ${silo_inc_suffixes})

    endif()
  endif()

  if (NOT Silo_INCLUDE_DIR)
    message(SEND_ERROR "Can not locate Silo include directory")
  endif()

  # Search for libraries 
  # Search order preference:
  #  (1) Silo_LIBRARY_DIR - check existence of path AND if the library file exists
  #  (2) Silo_DIR/<lib,Lib>
  #  (3) Default CMake paths See cmake --html-help=out.html file for more information.
  #
  set(silo_lib_names "siloh5")
  if (Silo_LIBRARY_DIR)
    if (EXISTS "${Silo_LIBRARY_DIR}")

      find_library(_Silo_LIBRARY
                   NAMES ${silo_lib_names}
                   HINTS ${Silo_LIBRARY_DIR}
                   NO_DEFAULT_PATH)
    else()
      message(SEND_ERROR "Silo_LIBRARY_DIR=${Silo_LIBRARY_DIR} does not exist")
      set(Silo_LIBRARY "Silo_LIBRARY-NOTFOUND")
    endif()

  else() 

    list(APPEND silo_lib_suffixes "lib" "Lib")
    if (Silo_DIR)
      if (EXISTS "${Silo_DIR}")

        find_library(_Silo_LIBRARY
                     NAMES ${silo_lib_names}
                     HINTS ${Silo_DIR}
                     PATH_SUFFIXES ${silo_lib_suffixes}
                     NO_DEFAULT_PATH)

      else()
        message(SEND_ERROR "Silo_DIR=${Silo_DIR} does not exist")
        set(_Silo_LIBRARY _Silo_LIBRARY-NOTFOUND)
      endif()    

    else()

      find_library(_Silo_LIBRARY
                   NAMES ${silo_lib_names}
                   PATH_SUFFIXES ${silo_lib_suffixes})

    endif()
  endif()

  # Create the target
  if (_Silo_LIBRARY)
    set(Silo_LIBRARY silo)
    add_imported_library(${Silo_LIBRARY}
                         LOCATION ${_Silo_LIBRARY}
                         LINK_LANGUAGES "C")
  else()			   
    message(SEND_ERROR "Can not locate Silo library")
  endif()    

  # Update the INCLUDE_DIRS and LIBRARIES variables
  set(Silo_INCLUDE_DIRS ${Silo_INCLUDE_DIR})
  set(Silo_LIBRARIES    ${Silo_LIBRARY})

  # Define the dependent libs
  set(_Silo_DEP_LIBS)

  # Silo depends on HDF5
  find_package(HDF5 QUIET REQUIRED)
  list(APPEND Silo_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS})			
  list(APPEND _Silo_DEP_LIBS ${HDF5_LIBRARY})

  set_target_properties(${Silo_LIBRARY} PROPERTIES
                        IMPORTED_LINK_INTERFACE_LIBRARIES "${_Silo_DEP_LIBS}")

endif(Silo_LIBRARIES AND Silo_INCLUDE_DIRS)  


# find_package_handle_standard_args should set Silo_FOUND but it does not!
if (Silo_LIBRARIES AND Silo_INCLUDE_DIRS)
  set(Silo_FOUND TRUE)
else()
  set(Silo_FOUND FALSE)
endif()

# Define the version
mark_as_advanced(
  Silo_INCLUDE_DIR
  Silo_INCLUDE_DIRS
  Silo_LIBRARY
  Silo_LIBRARIES
  Silo_LIBRARY_DIR
  Silo_UTILITIES
)
