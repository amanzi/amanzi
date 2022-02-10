# -*- mode: cmake -*-

#
# Amanzi SEACAS Find Module
#
# Usage:
#    Control the search through SEACAS_DIR or setting environment variable
#    SEACAS_ROOT to the SEACAS installation prefix.
#
#    This module does not search default paths! 
#
#    Following variables are set:
#    SEACAS_FOUND            (BOOL)       Flag indicating if SEACAS was found
#    SEACAS_INCLUDE_DIR      (PATH)       Path to the SEACAS include file
#    SEACAS_INCLUDE_DIRS     (LIST)       List of all required include files
#    SEACAS_LIBRARY_DIR      (PATH)       Path to the SEACAS library
#    SEACAS_LIBRARY          (FILE)       SEACAS library
#    SEACAS_LIBRARIES        (LIST)       List of all required SEACAS libraries
#
# #############################################################################

# Standard CMake modules see CMAKE_ROOT/Modules
include(FindPackageHandleStandardArgs)

include(AddImportedLibrary)

# Amanzi CMake functions see <root>/tools/cmake for source
include(PrintVariable)
include(AddPackageDependency)

if (SEACAS_LIBRARIES AND SEACAS_INCLUDE_DIRS)

  # Do nothing. Variables are set. No need to search again
  set(SEACAS_INCLUDE_DIR ${SEACAS_INCLUDE_DIRS})
  set(SEACAS_LIBRARY    ${SEACAS_LIBRARIES})

else()

  # Cache variables
  if (SEACAS_DIR)
    set(SEACAS_DIR "${SEACAS_DIR}" CACHE PATH "Path to search for SEACAS include and library files")
  endif()

  if (SEACAS_INCLUDE_DIR)
    set(SEACAS_INCLUDE_DIR "${SEACAS_INCLUDE_DIR}" CACHE PATH "Path to search for SEACAS include files")
  endif()

  if (SEACAS_LIBRARY_DIR)
    set(SEACAS_LIBRARY_DIR "${SEACAS_LIBRARY_DIR}" CACHE PATH "Path to search for SEACAS library files")
  endif()

  # Search for ExodusII include files
  set(exodus_inc_names "SEACASIoss_config.h")
  if (SEACAS_INCLUDE_DIR)
    if (EXISTS "${SEACAS_INCLUDE_DIR}")

      find_path(_exodusII_include_path
                NAMES ${exodus_inc_names}
                HINTS ${SEACAS_INCLUDE_DIR}
                NO_DEFAULT_PATH)

      if (NOT _exodusII_include_path)
        message(SEND_ERROR "Can not locate ${exodus_inc_names} in ${SEACAS_INCLUDE_DIR}")
      endif()

      set(SEACAS_INCLUDE_DIR "${_exodusII_include_path}")
    else()
      message(SEND_ERROR "SEACAS_INCLUDE_DIR=${SEACAS_INCLUDE_DIR} does not exist")
      set(SEACAS_INCLUDE_DIR "SEACAS_INCLUDE_DIR-NOTFOUND")
    endif()

  else() 

    set(exodus_inc_suffixes "include")
    if (SEACAS_DIR)
      if (EXISTS "${SEACAS_DIR}")

        find_path(SEACAS_INCLUDE_DIR
                  NAMES ${exodus_inc_names}
                  HINTS ${SEACAS_DIR}
                  PATH_SUFFIXES ${exodus_inc_suffixes}
                  NO_DEFAULT_PATH)

      else()
        message(SEND_ERROR "SEACAS_DIR=${SEACAS_DIR} does not exist")
        set(SEACAS_INCLUDE_DIR "SEACAS_INCLUDE_DIR-NOTFOUND")
      endif()    

    else()

      find_path(SEACAS_INCLUDE_DIR
                NAMES ${exodus_inc_names}
                PATH_SUFFIXES ${exodus_inc_suffixes})

    endif()
  endif()

  if (NOT SEACAS_INCLUDE_DIR )
    message(SEND_ERROR "Can not locate SEACAS include directory")
  endif()

  # Search for libraries 
  # set(exodus_lib_names "exoIIv2c" "exodus")
  set(exodus_lib_names "exodus")
  if (SEACAS_LIBRARY_DIR)
    if (EXISTS "${SEACAS_LIBRARY_DIR}")

      find_library(_SEACAS_LIBRARY
                   NAMES ${exodus_lib_names}
                   HINTS ${SEACAS_LIBRARY_DIR}
                   NO_DEFAULT_PATH)

      find_library(_SEACAS_Fortran_LIBRARY
                   NAMES exodus_for
                   HINTS ${SEACAS_LIBRARY_DIR}
                   NO_DEFAULT_PATH)

    else()
      message(SEND_ERROR "SEACAS_LIBRARY_DIR=${SEACAS_LIBRARY_DIR} does not exist")
      set(_SEACAS_LIBRARY "SEACAS_LIBRARY-NOTFOUND")
      set(_SEACAS_Fortran_LIBRARY "SEACAS_Fortran_LIBRARY-NOTFOUND")
    endif()

  else() 
    list(APPEND exodus_lib_suffixes "lib" "Lib")
    if (SEACAS_DIR)
      if (EXISTS "${SEACAS_DIR}")

        find_library(_SEACAS_LIBRARY
                     NAMES ${exodus_lib_names}
                     HINTS ${SEACAS_DIR}
                     PATH_SUFFIXES ${exodus_lib_suffixes}
                     NO_DEFAULT_PATH)
                
        find_library(_SEACAS_Fortran_LIBRARY
                     NAMES exodus_for
                     HINTS ${SEACAS_DIR}
                     PATH_SUFFIXES ${exodus_lib_suffixes}
                     NO_DEFAULT_PATH)

      else()
        message(SEND_ERROR "SEACAS_DIR=${SEACAS_DIR} does not exist")
        set(SEACAS_LIBRARY "SEACAS_LIBRARY-NOTFOUND")
        set(SEACAS_Fortran_LIBRARY "SEACAS_Fortran_LIBRARY-NOTFOUND")
      endif()    

    else()

      find_library(_SEACAS_LIBRARY
                   NAMES ${exodus_lib_names}
                   PATH_SUFFIXES ${exodus_lib_suffixes})

      find_library(_SEACAS_Fortran_LIBRARY
                   NAMES exodus_for
                   PATH_SUFFIXES ${exodus_lib_suffixes})
    endif()
  endif()

  # Create the library target store the name in SEACAS_LIBRARY
  if (_SEACAS_LIBRARY)
    set(SEACAS_LIBRARY exodus)
    add_imported_library(${SEACAS_LIBRARY}
                         LOCATION ${_SEACAS_LIBRARY}
                         LINK_LANGUAGES "C;CXX")
  else()
    message(SEND_ERROR "Can not locate SEACAS library")
  endif()

  if (_SEACAS_Fortran_LIBRARY)
    set(SEACAS_Fortran_LIBRARY exodusii_for)
    add_imported_library(${SEACAS_Fortran_LIBRARY}
                         LOCATION ${_SEACAS_Fortran_LIBRARY}
                         LINK_LANGUAGES "Fortran")
    set_target_properties(${SEACAS_Fortran_LIBRARY} PROPERTIES
                          INTERFACE_LINK_LIBRARIES "${Exodus_LIBRARY}")
  endif()

  # Define prerequisite packages
  set(SEACAS_INCLUDE_DIRS ${SEACAS_INCLUDE_DIR})
  set(SEACAS_LIBRARIES    ${SEACAS_LIBRARY})

  # Search for NetCDF
  find_package(NetCDF QUIET REQUIRED)
  set_target_properties(${SEACAS_LIBRARY} PROPERTIES
                        INTERFACE_LINK_LIBRARIES "${NetCDF_C_LIBRARIES}")
  list(APPEND SEACAS_INCLUDE_DIRS ${NetCDF_INCLUDE_DIRS})

endif(SEACAS_LIBRARIES AND SEACAS_INCLUDE_DIRS )    

# Send useful message if everything is found
find_package_handle_standard_args(SEACAS DEFAULT_MSG
                                  SEACAS_INCLUDE_DIR
                                  SEACAS_LIBRARIES)

# find_package should set SEACAS_FOUND but it does not!
if (SEACAS_LIBRARIES AND SEACAS_INCLUDE_DIRS)
  set(SEACAS_FOUND TRUE)
else()
  set(SEACAS_FOUND FALSE)
endif()

mark_as_advanced(
  SEACAS_VERSION
  SEACAS_INCLUDE_DIR
  SEACAS_INCLUDE_DIRS
  SEACAS_LIBRARY
  SEACAS_LIBRARIES
  SEACAS_LIBRARY_DIR
)
