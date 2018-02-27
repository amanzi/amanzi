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

if (ALQUIMIA_LIBRARIES AND ALQUIMIA_INCLUDE_DIRS)

  # Do nothing. Variables are set. No need to search again

elseif(ALQUIMIA_DIR)

  set(ALQUIMIA_INCLUDE_DIR ${ALQUIMIA_DIR}/include)
  set(ALQUIMIA_LIBRARY_DIR ${ALQUIMIA_DIR}/lib)
  set(ALQUIMIA_TARGET alquimia)

  find_library(_ALQUIMIA_LIBRARY
               NAMES alquimia
               PATHS ${ALQUIMIA_LIBRARY_DIR})

  if (_ALQUIMIA_LIBRARY)
    add_imported_library(${ALQUIMIA_TARGET} 
                         LOCATION ${_ALQUIMIA_LIBRARY}
                         LINK_LANGUAGES "C")
    set(ALQUIMIA_LIBRARY ${ALQUIMIA_TARGET})
  endif()    

  # Define the LIBRARIES and INCLUDE_DIRS
  set(ALQUIMIA_INCLUDE_DIRS ${ALQUIMIA_INCLUDE_DIR})
  set(ALQUIMIA_LIBRARIES ${ALQUIMIA_LIBRARY} ${CRUNCHTOPE_LIBRARY} ${PFLOTRAN_LIBRARIES})

endif(ALQUIMIA_LIBRARIES AND ALQUIMIA_INCLUDE_DIRS)    

# Send useful message if everything is found
find_package_handle_standard_args(ALQUIMIA DEFAULT_MSG
                                  ALQUIMIA_INCLUDE_DIRS
                                  ALQUIMIA_LIBRARIES)

# find_package)handle)standard_args should set ALQUIMIA_FOUND but it does not!
if (ALQUIMIA_LIBRARIES AND ALQUIMIA_INCLUDE_DIRS)
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
