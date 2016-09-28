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

else()
    set(ALQUIMIA_INCLUDE_DIR ${TPL_INSTALL_PREFIX}/include/alquimia)
    set(ALQUIMIA_LIBRARY_DIR ${TPL_INSTALL_PREFIX}/lib)
    set(ALQUIMIA_LIBRARY libalquimia.a)
    find_library(_ALQUIMIA_CRUNCH_LIBRARY
                 NAMES "libcrunchchem.a"
                 PATHS $ENV{CRUNCH_DIR})

    if ( _ALQUIMIA_CRUNCH_LIBRARY )
        add_imported_library(${ALQUIMIA_CRUNCH_TARGET}
	                     LOCATION ${_ALQUIMIA_CRUNCH_LIBRARY}
                             LINK_LANGUAGES "Fortran")
      set(ALQUIMIA_CRUNCH_LIBRARY ${ALQUIMIA_CRUNCH_TARGET})
    endif()    
    
    # Define the LIBRARIES and INCLUDE_DIRS
    set(ALQUIMIA_INCLUDE_DIRS ${ALQUIMIA_INCLUDE_DIR})
    set(ALQUIMIA_LIBRARIES ${ALQUIMIA_LIBRARY} ${ALQUIMIA_CRUNCH_LIBRARY} ${PFLOTRAN_LIBRARIES})
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
  ALQUIMIA_LIBRARY_DIR
  ALQUIMIA_LIBRARY
  ALQUIMIA_LIBRARIES
)
