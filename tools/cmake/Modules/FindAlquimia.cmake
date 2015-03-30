# -*- mode: cmake -*-

#
# Amanzi ALQUIMIA Find Module
#
# Usage:
#    Control the search through ALQUIMIA_DIR or setting environment variable
#    ALQUIMIA_ROOT to the ALQUIMIA installation prefix.
#
#    This module does not search default paths! 
#
#    Following variables are set:
#    ALQUIMIA_FOUND            (BOOL)       Flag indicating if ALQUIMIA was found
#    ALQUIMIA_INCLUDE_DIR      (PATH)       Path to the ALQUIMIA include file
#    ALQUIMIA_INCLUDE_DIRS     (LIST)       List of all required include files
#    ALQUIMIA_LIBRARY_DIR      (PATH)       Path to the ALQUIMIA library
#    ALQUIMIA_LIBRARIES        (LIST)       List of all required ALQUIMIA libraries
#
# #############################################################################

# Standard CMake modules see CMAKE_ROOT/Modules
include(FindPackageHandleStandardArgs)

# Amanzi CMake functions see <root>/tools/cmake for source
include(PrintVariable)
include(AddPackageDependency)

if ( ALQUIMIA_LIBRARIES AND ALQUIMIA_INCLUDE_DIRS )

    # Do nothing. Variables are set. No need to search again

else(ALQUIMIA_LIBRARIES AND ALQUIMIA_INCLUDE_DIRS)

    message(STATUS "ALQUIMIA_DIR = ${ALQUIMIA_DIR}")
    # Cache variables
    if(ALQUIMIA_DIR)
      set(ALQUIMIA_DIR "${ALQUIMIA_DIR}" CACHE PATH "Path to search for ALQUIMIA include and library files")
    endif()

    if(ALQUIMIA_INCLUDE_DIR)
      set(ALQUIMIA_INCLUDE_DIR ${ALQUIMIA_INCLUDE_DIR} CACHE PATH "Path to search for ALQUIMIA include files")
    else()
      find_path(ALQUIMIA_INCLUDE_DIR alquimia_interface.h ${ALQUIMIA_DIR}/include)
      if ( NOT ALQUIMIA_INCLUDE_DIR )
        message(SEND_ERROR "Cannot locate ALQUIMIA include directory")
      endif()
    endif()

    if(ALQUIMIA_LIBRARY_DIR)
      set(ALQUIMIA_LIBRARY_DIR "${ALQUIMIA_LIBRARY_DIR}" CACHE PATH "Path to search for ALQUIMIA library files")
    else()
      find_path(ALQUIMIA_LIBRARY_DIR NAMES libalquimia_c.a PATHS ${ALQUIMIA_DIR}/lib)
      if ( NOT ALQUIMIA_LIBRARY_DIR )
        message(SEND_ERROR "Cannot locate ALQUIMIA library directory")
      endif()
    endif()

    # Search for libraries 

    set( ALQUIMIA_C_TARGET alquimia_c )
    set( ALQUIMIA_CUTILS_TARGET alquimia_cutils )
    set( ALQUIMIA_FORTRAN_TARGET alquimia_fortran )
    set( ALQUIMIA_CRUNCH_TARGET "libcrunchchem.a" )

    find_library(_ALQUIMIA_C_LIBRARY
                 NAMES alquimia_c
                 PATHS ${ALQUIMIA_LIBRARY_DIR})

    find_library(_ALQUIMIA_CUTILS_LIBRARY
                 NAMES alquimia_cutils
                 PATHS ${ALQUIMIA_LIBRARY_DIR})

    find_library(_ALQUIMIA_FORTRAN_LIBRARY
                 NAMES alquimia_fortran
                 PATHS ${ALQUIMIA_LIBRARY_DIR})

    find_library(_ALQUIMIA_CRUNCH_LIBRARY
                 NAMES "libcrunchchem.a"
                 PATHS $ENV{CRUNCH_DIR})

    if ( _ALQUIMIA_C_LIBRARY )
      add_imported_library(${ALQUIMIA_C_TARGET} 
                           LOCATION ${_ALQUIMIA_C_LIBRARY}
                           LINK_LANGUAGES "C")
      set(ALQUIMIA_C_LIBRARY ${ALQUIMIA_C_TARGET})
    endif()    

    if ( _ALQUIMIA_CUTILS_LIBRARY )
        add_imported_library(${ALQUIMIA_CUTILS_TARGET}
	                     LOCATION ${_ALQUIMIA_CUTILS_LIBRARY}
                             LINK_LANGUAGES "C")
      set(ALQUIMIA_CUTILS_LIBRARY ${ALQUIMIA_CUTILS_TARGET})
    endif()    

    if ( _ALQUIMIA_FORTRAN_LIBRARY )
        add_imported_library(${ALQUIMIA_FORTRAN_TARGET}
	                     LOCATION ${_ALQUIMIA_FORTRAN_LIBRARY}
                             LINK_LANGUAGES "Fortran")
      set(ALQUIMIA_FORTRAN_LIBRARY ${ALQUIMIA_FORTRAN_TARGET})
    endif()

    if ( _ALQUIMIA_CRUNCH_LIBRARY )
        add_imported_library(${ALQUIMIA_CRUNCH_TARGET}
	                     LOCATION ${_ALQUIMIA_CRUNCH_LIBRARY}
                             LINK_LANGUAGES "Fortran")
      set(ALQUIMIA_CRUNCH_LIBRARY ${ALQUIMIA_CRUNCH_TARGET})
    endif()    
    
    # Define the LIBRARIES and INCLUDE_DIRS
    set(ALQUIMIA_INCLUDE_DIRS ${ALQUIMIA_INCLUDE_DIR})
    set(ALQUIMIA_LIBRARIES    
      ${ALQUIMIA_C_LIBRARY} 
      ${ALQUIMIA_CUTILS_LIBRARY} 
      ${ALQUIMIA_FORTRAN_LIBRARY}
      ${ALQUIMIA_CRUNCH_LIBRARY} )

endif(ALQUIMIA_LIBRARIES AND ALQUIMIA_INCLUDE_DIRS )    

# Send useful message if everything is found
find_package_handle_standard_args(ALQUIMIA DEFAULT_MSG
                                           ALQUIMIA_INCLUDE_DIRS
					   ALQUIMIA_LIBRARIES)

# find_package)handle)standard_args should set ALQUIMIA_FOUND but it does not!
if ( ALQUIMIA_LIBRARIES AND ALQUIMIA_INCLUDE_DIRS)
    set(ALQUIMIA_FOUND TRUE)
else()
    set(ALQUIMIA_FOUND FALSE)
endif()

mark_as_advanced(
  ALQUIMIA_INCLUDE_DIR
  ALQUIMIA_INCLUDE_DIRS
  ALQUIMIA_LIBRARY
  ALQUIMIA_LIBRARIES
  ALQUIMIA_LIBRARY_DIR
)
