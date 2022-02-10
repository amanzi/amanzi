# -*- mode: cmake -*-

#
# Amanzi XERCES Find Module
#
# Usage:
#    Control the search through XERCES_DIR or setting environment variable
#    XERCES_ROOT to the XERCES installation prefix.
#
#    This module does not search default paths! 
#
#    Following variables are set:
#    XERCES_FOUND            (BOOL)       Flag indicating if XERCES was found
#    XERCES_INCLUDE_DIR      (PATH)       Path to the XERCES include file
#    XERCES_INCLUDE_DIRS     (LIST)       List of all required include files
#    XERCES_LIBRARY_DIR      (PATH)       Path to the XERCES library
#    XERCES_LIBRARY          (FILE)       XERCES library
#    XERCES_LIBRARIES        (LIST)       List of all required XERCES libraries
#
# #############################################################################

# Standard CMake modules see CMAKE_ROOT/Modules
include(FindPackageHandleStandardArgs)

# Amanzi CMake functions see <root>/tools/cmake for source
include(PrintVariable)
include(AddPackageDependency)

if (XERCES_LIBRARIES AND XERCES_INCLUDE_DIRS)

  # Do nothing. Variables are set. No need to search again

else(XERCES_LIBRARY_DIR AND XERCES_INCLUDE_DIR)

  # Cache variables
  if (XERCES_DIR)
    set(XERCES_DIR "${XERCES_DIR}" CACHE PATH "Path to search for XERCES include and library files")
  endif()

  if (XERCES_INCLUDE_DIR)
    set(XERCES_INCLUDE_DIR "${XERCES_INCLUDE_DIR}" CACHE PATH "Path to search for XERCES include files")
  endif()

  if (XERCES_LIBRARY_DIR)
    set(XERCES_LIBRARY_DIR "${XERCES_LIBRARY_DIR}" CACHE PATH "Path to search for XERCES library files")
  endif()

  # Search for include files
  # Search order preference:
  #  (1) XERCES_INCLUDE_DIR - check existence of path AND if the include files exist
  #  (2) XERCES_DIR/<include>
  #  (3) Default CMake paths See cmake --html-help=out.html file for more information.
  #
  set(xerces_inc_names "dom")
  #set(xerces_inc_suffixes  "xercesc" )
  if (XERCES_INCLUDE_DIR)
    set(DOM_INCLUDE_DIR "${XERCES_INCLUDE_DIR}/dom")
    if (EXISTS "${DOM_INCLUDE_DIR}")
      set(XERCES_INCLUDE_DIR "${XERCES_INCLUDE_DIR}")
    else()
      set(XERCES_INCLUDE_DIR "XERCES_INCLUDE_DIR-NOTFOUND")
    endif()
  else() 

  set(xerces_inc_suffixes "include" "include/xercesc" )
    if (XERCES_DIR)
      if (EXISTS "${XERCES_DIR}")
        message (STATUS "EIB >>> trying to find ${xerces_inc_names} in "
                        " ${XERCES_DIR} with ${xerces_inc_suffixes}")

        find_path(XERCES_INCLUDE_DIR
                  NAMES ${xerces_inc_names}
	          HINTS ${XERCES_DIR}
                  PATH_SUFFIXES ${xerces_inc_suffixes}
	          NO_DEFAULT_PATH)

      else()
        message(SEND_ERROR "XERCES_DIR=${XERCES_DIR} does not exist")
        set(XERCES_INCLUDE_DIR "XERCES_INCLUDE_DIR-NOTFOUND")
      endif()    

    else()
      message(STATUS "EIB >>> last option look in ${xerces_inc_suffixes} "
                     "for ${xerces_inc_names}")
      find_path(XERCES_INCLUDE_DIR
                NAMES ${xerces_inc_names}
                PATH_SUFFIXES ${xerces_inc_suffixes})
    endif()
  endif()

  if (NOT XERCES_INCLUDE_DIR)
    message(SEND_ERROR "Can not locate XERCES include directory")
  endif()

  # Search for libraries 
  # Search order preference:
  #  (1) XERCES_LIBRARY_DIR - check existence of path AND if the include files exist
  #  (2) XERCES_DIR/<lib,Lib>
  #  (3) Default CMake paths See cmake --html-help=out.html file for more information.
  #
  set(xerces_lib_names "xerces-c")
  set(xerces_lib_suffixes "xerces-c")
  if (XERCES_LIBRARY_DIR)
    if (EXISTS "${XERCES_LIBRARY_DIR}")

      find_library(XERCES_LIBRARY
                   NAMES ${xerces_lib_names}
                   HINTS ${XERCES_LIBRARY_DIR}
                   PATH_SUFFIXES ${xerces_lib_suffixes}
                   NO_DEFAULT_PATH)
    else()
      message(SEND_ERROR "XERCES_LIBRARY_DIR=${XERCES_LIBRARY_DIR} does not exist")
      set(XERCES_LIBRARY "XERCES_LIBRARY-NOTFOUND")
    endif()

  else() 
    list(APPEND xerces_lib_suffixes "lib" "lib/xerces")
    list(APPEND xerces_Lib_suffixes "Lib" "Lib/xerces")
    if (XERCES_DIR)
      if (EXISTS "${XERCES_DIR}")
	message (STATUS "EIB >>> trying to find ${xerces_lib_names} in "
		        " ${XERCES_DIR} with ${xerces_lib_suffixes}")

        find_library(XERCES_LIBRARY
                     NAMES ${xerces_lib_names}
                     HINTS ${XERCES_DIR}
                     PATH_SUFFIXES ${xerces_lib_suffixes}
                     NO_DEFAULT_PATH)

      else()
        message(SEND_ERROR "XERCES_DIR=${XERCES_DIR} does not exist")
        set(XERCES_LIBRARY "XERCES_LIBRARY-NOTFOUND")
      endif()    

    else()

      find_library(XERCES_LIBRARY
                   NAMES ${xerces_lib_names}
                   PATH_SUFFIXES ${xerces_lib_suffixes})

    endif()
  endif()

  if (NOT XERCES_LIBRARY )
    set(XERCES_LIBRARY XERCES_LIBRARY-NOTFOUND)
    message(SEND_ERROR "Can not locate XERCES library")
  endif()    

  # Grab the library dependencies from the setting file xerces-c.pc
  #file (STRINGS ${XERCES_DIR}/lib/pkgconfig/xerces-c.pc XERCES_LIBRARY_DEPS REGEX "Libs.private:*")

  if (NOT "${XERCES_LIBRARY_DEPS}" STREQUAL "")
    STRING(REGEX REPLACE "Libs.private:" "" XERCES_LIBRARY_DEPS ${XERCES_LIBRARY_DEPS})
    STRING(REGEX REPLACE " " "" XERCES_LIBRARY_DEPS "${XERCES_LIBRARY_DEPS}")
  endif()

  # For now we don't recurse on *.la files
  set(XERCES_ICU_LIBRARIES "")
  foreach(ln ${XERCES_LIBRARY_DEPS})
    STRING(REGEX MATCH "\\.la" OUT_libtool ${ln})
    # Drop system libraries (-L/usr/*) because they should be in the system linker already
    # -- would be more robust to get ld search path, and then drop overlapping.
    STRING(REGEX MATCH "[-][L]/usr" OUT_lib_system ${ln})
    if (NOT OUT_libtool AND NOT OUT_lib_system) 
      list(APPEND XERCES_ICU_LIBRARIES ${ln})
    endif()
  endforeach()

  # Add dependency on frameworks for OSX (it doesn't appear in this list for some reason)
  if (APPLE) 
    list (APPEND XERCES_ICU_LIBRARIES "-framework CoreServices")
  endif()

  # Define the LIBRARIES and INCLUDE_DORS
  set(XERCES_INCLUDE_DIRS ${XERCES_INCLUDE_DIR})
  set(XERCES_LIBRARIES ${XERCES_LIBRARY})
  list(APPEND XERCES_LIBRARIES ${XERCES_LIBRARY_DEPS})

endif(XERCES_LIBRARIES AND XERCES_INCLUDE_DIRS )    

# Send useful message if everything is found
find_package_handle_standard_args(XERCES DEFAULT_MSG
                                  XERCES_INCLUDE_DIRS
                                  XERCES_LIBRARIES)

# find_package)handle)standard_args should set XERCES_FOUND but it does not!
if (XERCES_LIBRARIES AND XERCES_INCLUDE_DIRS)
  set(XERCES_FOUND TRUE)
else()
  set(XERCES_FOUND FALSE)
endif()

mark_as_advanced(
  XERCES_INCLUDE_DIR
  XERCES_INCLUDE_DIRS
  XERCES_LIBRARY
  XERCES_LIBRARIES
  XERCES_LIBRARY_DIR
  XERCES_ICU_LIBRARIES
)
