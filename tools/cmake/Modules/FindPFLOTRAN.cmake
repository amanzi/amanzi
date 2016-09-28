# -*- mode: cmake -*-

#
# Amanzi PFLOTRAN Find Module
#
# Usage:
#    Control the search through PFLOTRAN_DIR or setting environment variable
#    PFLOTRAN_ROOT to the PFLOTRAN installation prefix.
#
#    This module does not search default paths! 
#
#    Following variables are set:
#    PFLOTRAN_FOUND            (BOOL)       Flag indicating if PFLOTRAN was found
#    PFLOTRAN_INCLUDE_DIR      (PATH)       Path to the PFLOTRAN include file
#    PFLOTRAN_INCLUDE_DIRS     (LIST)       List of all required include files
#    PFLOTRAN_LIBRARY_DIR      (PATH)       Path to the PFLOTRAN library
#    PFLOTRAN_LIBRARY          (FILE)       PFLOTRAN library
#    PFLOTRAN_LIBRARIES        (LIST)       List of all required PFLOTRAN libraries
#
# #############################################################################

# Standard CMake modules see CMAKE_ROOT/Modules
include(FindPackageHandleStandardArgs)

# Amanzi CMake functions see <root>/tools/cmake for source
include(PrintVariable)
include(AddPackageDependency)

if (PFLOTRAN_LIBRARIES)

    # Do nothing. Variables are set. No need to search again

else()

    # Cache variables
    if(PFLOTRAN_DIR)
      set(PFLOTRAN_DIR "${PFLOTRAN_DIR}" CACHE PATH "Path to search for PFLOTRAN include and library files")
    endif()

    if(PFLOTRAN_LIBRARY_DIR)
        set(PFLOTRAN_LIBRARY_DIR "${PFLOTRAN_LIBRARY_DIR}" CACHE PATH "Path to search for PFLOTRAN library files")
    else()
      find_path(PFLOTRAN_LIBRARY_DIR NAMES libpflotranchem.a PATHS ${PFLOTRAN_DIR}/src/pflotran)
      if ( NOT PFLOTRAN_LIBRARY_DIR )
        message(SEND_ERROR "Cannot locate PFLTORAN library directory")
      endif()
    endif()

    # Search for libraries 
    
    set( PFLOTRAN_TARGET pflotranchem )

    find_library(_PFLOTRAN_LIBRARY
                 NAMES pflotranchem
                 PATHS ${PFLOTRAN_LIBRARY_DIR}
                 NO_DEFAULT_PATH )

    if ( _PFLOTRAN_LIBRARY ) 
      add_imported_library(${PFLOTRAN_TARGET}
	                   LOCATION ${_PFLOTRAN_LIBRARY}
                           LINK_LANGUAGES "Fortran")
      set(PFLOTRAN_LIBRARY ${PFLOTRAN_TARGET})
    else()
      set(PFLOTRAN_LIBRARY PFLOTRAN_LIBRARY-NOTFOUND)
      message(SEND_ERROR "Cannot locate PFLOTRAN library")
    endif()    
    
   
    # Define the LIBRARIES and INCLUDE_DIRS
    set(PFLOTRAN_LIBRARIES    ${PFLOTRAN_LIBRARY})

endif()    

# Send useful message if everything is found
find_package_handle_standard_args(PFLOTRAN DEFAULT_MSG
                                           PFLOTRAN_INCLUDE_DIRS
					   PFLOTRAN_LIBRARIES)

# find_package)handle)standard_args should set PFLOTRAN_FOUND but it does not!
if (PFLOTRAN_LIBRARIES)
    set(PFLOTRAN_FOUND TRUE)
else()
    set(PFLOTRAN_FOUND FALSE)
endif()

mark_as_advanced(
  PFLOTRAN_LIBRARIES
  PFLOTRAN_LIBRARY_DIR
)
