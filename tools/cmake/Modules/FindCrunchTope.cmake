# -*- mode: cmake -*-

#
# Amanzi CRUNCHTOPE Find Module
#
# Usage:
#    Control the search through CRUNCHTOPE_DIR installation prefix.
#
#    This module does not search default paths! 
#
#    Following variables are set:
#    CRUNCHTOPE_FOUND          (BOOL)    Flag indicating if CRUNCHTOPE was found
#    CRUNCHTOPE_INCLUDE_DIR    (PATH)    Path to the CRUNCHTOPE include file
#    CRUNCHTOPE_INCLUDE_DIRS   (LIST)    List of all required include files
#    CRUNCHTOPE_LIBRARY_DIR    (PATH)    Path to the CRUNCHTOPE library
#    CRUNCHTOPE_LIBRARIES      (LIST)    List of all required CRUNCHTOPE libraries
#
# #############################################################################

# Standard CMake modules see CMAKE_ROOT/Modules
include(FindPackageHandleStandardArgs)

# Amanzi CMake functions see <root>/tools/cmake for source
include(PrintVariable)
include(AddPackageDependency)

if ( CRUNCHTOPE_LIBRARIES AND CRUNCHTOPE_INCLUDE_DIRS )

    # Do nothing. Variables are set. No need to search again

elseif ( CRUNCHTOPE_DIR )

    set(CRUNCHTOPE_INCLUDE_DIR ${CRUNCHTOPE_DIR}/lib)
    set(CRUNCHTOPE_LIBRARY_DIR ${CRUNCHTOPE_DIR}/lib)
    set(CRUNCHTOPE_TARGET crunchchem.a)

    find_library(_CRUNCHTOPE_LIBRARY
                 NAMES crunchchem
                 PATHS ${CRUNCHTOPE_DIR}/lib)

    if ( _CRUNCHTOPE_LIBRARY )
        add_imported_library(${CRUNCHTOPE_TARGET}
	                     LOCATION ${_CRUNCHTOPE_LIBRARY}
                             LINK_LANGUAGES "Fortran")
        set(CRUNCHTOPE_LIBRARY ${CRUNCHTOPE_TARGET})
    endif()    
    
    # Define the LIBRARIES and INCLUDE_DIRS
    set(CRUNCHTOPE_INCLUDE_DIRS ${CRUNCHTOPE_INCLUDE_DIR})
    set(CRUNCHTOPE_LIBRARIES ${CRUNCHTOPE_LIBRARY})

endif( CRUNCHTOPE_LIBRARIES AND CRUNCHTOPE_INCLUDE_DIRS )

# Send useful message if everything is found
find_package_handle_standard_args(CRUNCHTOPE DEFAULT_MSG
                                             CRUNCHTOPE_INCLUDE_DIRS
					     CRUNCHTOPE_LIBRARIES)

# find_package should set CRUNCHTOPE_FOUND but it does not!
if (CRUNCHTOPE_LIBRARIES AND CRUNCHTOPE_INCLUDE_DIRS)
    set(CRUNCHTOPE_FOUND TRUE)
else()
    set(CRUNCHTOPE_FOUND FALSE)
endif()

mark_as_advanced(
  CRUNCHTOPE_INCLUDE_DIR
  CRUNCHTOPE_LIBRARY_DIR
  CRUNCHTOPE_LIBRARIES
)
