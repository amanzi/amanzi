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

# Standard CMake modules see CMAKE_ROOT/Modules
include(FindPackageHandleStandardArgs)

# Amanzi CMake functions see <root>/tools/cmake for source
include(PrintVariable)
include(AddPackageDependency)

if (Alquimia_LIBRARIES AND Alquimia_INCLUDE_DIRS)

  # Do nothing. Variables are set. No need to search again

elseif(Alquimia_DIR)

  set(Alquimia_INCLUDE_DIR ${Alquimia_DIR}/include)
  set(Alquimia_LIBRARY_DIR ${Alquimia_DIR}/lib)
  set(Alquimia_TARGET alquimia)

  find_library(_Alquimia_LIBRARY
               NAMES alquimia
               PATHS ${Alquimia_LIBRARY_DIR})

  if (_Alquimia_LIBRARY)
    add_imported_library(${Alquimia_TARGET} 
                         LOCATION ${_Alquimia_LIBRARY}
                         LINK_LANGUAGES "C")
    set(Alquimia_LIBRARY ${Alquimia_TARGET})
  endif()    

  # Define the LIBRARIES and INCLUDE_DIRS
  set(Alquimia_INCLUDE_DIRS ${Alquimia_INCLUDE_DIR})
  set(Alquimia_LIBRARIES ${Alquimia_LIBRARY} ${CrunchTope_LIBRARY} ${PFLOTRAN_LIBRARIES})

endif(Alquimia_LIBRARIES AND Alquimia_INCLUDE_DIRS)    

# Send useful message if everything is found
find_package_handle_standard_args(Alquimia DEFAULT_MSG
                                  Alquimia_INCLUDE_DIRS
                                  Alquimia_LIBRARIES)

# find_package)handle)standard_args should set Alquimia_FOUND but it does not!
if (Alquimia_LIBRARIES AND Alquimia_INCLUDE_DIRS)
  set(Alquimia_FOUND TRUE)
else()
  set(Alquimia_FOUND FALSE)
endif()

mark_as_advanced(
  Alquimia_INCLUDE_DIR
  Alquimia_LIBRARY_DIR
  Alquimia_LIBRARY
  Alquimia_LIBRARIES
)
