# -*- mode: cmake -*-

#
# Amanzi Ginkgo Find Module
#
# Usage:
#    Control the search through Ginkgo_DIR or setting environment variable
#    Ginkgo_ROOT to the Ginkgo installation prefix.
#
#    This module does not search default paths! 
#
#    Following variables are set:
#    Ginkgo_FOUND            (BOOL)       Flag indicating if Ginkgo was found
#    Ginkgo_INCLUDE_DIR      (PATH)       Path to the Ginkgo include file
#    Ginkgo_INCLUDE_DIRS     (LIST)       List of all required include files
#    Ginkgo_LIBRARY_DIR      (PATH)       Path to the Ginkgo library
#    Ginkgo_LIBRARY          (FILE)       Ginkgo library
#    Ginkgo_LIBRARIES        (LIST)       List of all required Ginkgo libraries
#
# #############################################################################

# Standard CMake modules see CMAKE_ROOT/Modules
#include(FindPackageHandleStandardArgs)

# Amanzi CMake functions see <root>/tools/cmake for source
include(PrintVariable)
include(AddPackageDependency)

if (Ginkgo_LIBRARIES AND Ginkgo_INCLUDE_DIRS)

  # Do nothing. Variables are set. No need to search again

else(Ginkgo_LIBRARIES AND Ginkgo_INCLUDE_DIRS)

  # Cache variables
  if (Ginkgo_DIR)
    set(Ginkgo_DIR "${Ginkgo_DIR}" CACHE PATH "Path to search for Ginkgo include and library files")
  endif()

  if (Ginkgo_INCLUDE_DIR)
    set(Ginkgo_INCLUDE_DIR "${Ginkgo_INCLUDE_DIR}" CACHE PATH "Path to search for Ginkgo include files")
  endif()

  if (Ginkgo_LIBRARY_DIR)
    set(Ginkgo_LIBRARY_DIR "${Ginkgo_LIBRARY_DIR}" CACHE PATH "Path to search for Ginkgo library files")
  endif()

  set(Ginkgo_inc_names "ginkgo/ginkgo.hpp")
  if (Ginkgo_INCLUDE_DIR)
    # Search for include files in provided Ginkgo_INCLUDE_DIR
    if (EXISTS "${Ginkgo_INCLUDE_DIR}")

      find_path(Ginkgo_include_path
                NAMES ${Ginkgo_inc_names}
                HINTS ${Ginkgo_INCLUDE_DIR}
                NO_DEFAULT_PATH)

      if (NOT Ginkgo_include_path)
        message(SEND_ERROR "Can not locate ${Ginkgo_inc_names} in ${Ginkgo_INCLUDE_DIR}")
      endif()

    else()
      message(SEND_ERROR "Ginkgo_INCLUDE_DIR=${Ginkgo_INCLUDE_DIR} does not exist")
      set(Ginkgo_INCLUDE_DIR "Ginkgo_INCLUDE_DIR-NOTFOUND")
    endif()

  else() 

    # Search for include files in Ginkgo_DIR/include
    set(Ginkgo_inc_suffixes "include/ginkgp")
    if (Ginkgo_DIR)
      if (EXISTS "${Ginkgo_DIR}")

        find_path(Ginkgo_INCLUDE_DIR
                  NAMES ${Ginkgo_inc_names}
                  HINTS ${Ginkgo_DIR}
                  PATH_SUFFIXES ${Ginkgo_inc_suffixes}
                  NO_DEFAULT_PATH)

      else()
        message(SEND_ERROR "Ginkgo_DIR=${Ginkgo_DIR} does not exist")
        set(Ginkgo_INCLUDE_DIR "Ginkgo_INCLUDE_DIR-NOTFOUND")
      endif()    

    else()

      # Search for include files in default path
      find_path(Ginkgo_INCLUDE_DIR
                NAMES ${Ginkgo_inc_names}
                PATH_SUFFIXES ${Ginkgo_inc_suffixes})

    endif()
  endif()

  if (NOT Ginkgo_INCLUDE_DIR)
    message(SEND_ERROR "Can not locate Ginkgo include directory")
  endif()

  # Search for libraries 
  set(Ginkgo_lib_names "ginkgo")
  if (Ginkgo_LIBRARY_DIR)
    # Search for library files in provided Ginkgo_LIBRARY_DIR
    if (EXISTS "${Ginkgo_LIBRARY_DIR}")

      find_library(_Ginkgo_LIBRARY
                   NAMES ${Ginkgo_lib_names}
                   HINTS ${Ginkgo_LIBRARY_DIR}
                   NO_DEFAULT_PATH)
    else()
      message(SEND_ERROR "Ginkgo_LIBRARY_DIR=${Ginkgo_LIBRARY_DIR} does not exist")
      set(Ginkgo_LIBRARY "Ginkgo_LIBRARY-NOTFOUND")
    endif()

  else() 

    # Search for library files in Ginkgo_DIR/lib64
    list(APPEND Ginkgo_lib_suffixes "lib" "Lib" "lib64")
    if (Ginkgo_DIR)
      if (EXISTS "${Ginkgo_DIR}" )

         find_library(_Ginkgo_LIBRARY
                      NAMES ${Ginkgo_lib_names}
                      HINTS ${Ginkgo_DIR}
                      PATH_SUFFIXES ${Ginkgo_lib_suffixes}
                      NO_DEFAULT_PATH)

      else()
        message(SEND_ERROR "Ginkgo_DIR=${Ginkgo_DIR} does not exist")
        set(_Ginkgo_LIBRARY _Ginkgo_LIBRARY-NOTFOUND)
      endif()    

    else()

      # Search for library files in default path
      find_library(_Ginkgo_LIBRARY
                   NAMES ${Ginkgo_lib_names}
                   PATH_SUFFIXES ${Ginkgo_lib_suffixes})

    endif()
  endif()

  # Create the target
  if (_Ginkgo_LIBRARY)
    set(Ginkgo_LIBRARY ${Ginkgo_lib_names})
    get_filename_component(Ginkgo_LIBRARY_DIR ${_Ginkgo_LIBRARY} DIRECTORY)

    add_imported_library(${Ginkgo_LIBRARY}
                         LOCATION ${_Ginkgo_LIBRARY}
                   LINK_LANGUAGES "C")
  else()         
    message(SEND_ERROR "Can not locate Ginkgo library")
  endif()    

  # Update the INCLUDE_DIRS and LIBRARIES variables
  set(Ginkgo_INCLUDE_DIRS ${Ginkgo_INCLUDE_DIR})
  set(Ginkgo_LIBRARIES    ${_Ginkgo_LIBRARY})


endif(Ginkgo_LIBRARIES AND Ginkgo_INCLUDE_DIRS)    

# Send useful message if everything is found
 find_package_handle_standard_args(Ginkgo DEFAULT_MSG
                                   Ginkgo_LIBRARIES
                                   Ginkgo_INCLUDE_DIRS)

set(Ginkgo_FOUND TRUE)

