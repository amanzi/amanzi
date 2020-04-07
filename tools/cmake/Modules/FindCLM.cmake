# -*- mode: cmake -*-

#
# Amanzi CLM Find Module
#
# Usage:
#    Control the search through CLM_DIR installation prefix.
#
#    This module does not search default paths! 
#
#    Following variables are set:
#    CLM_FOUND          (BOOL)    Flag indicating if CLM was found
#    CLM_INCLUDE_DIR    (PATH)    Path to the CLM include file
#    CLM_INCLUDE_DIRS   (LIST)    List of all required include files
#    CLM_LIBRARY_DIR    (PATH)    Path to the CLM library
#    CLM_LIBRARIES      (LIST)    List of all required CLM libraries
#
# #############################################################################

# Standard CMake modules see CMAKE_ROOT/Modules
include(FindPackageHandleStandardArgs)

# Amanzi CMake functions see <root>/tools/cmake for source
include(PrintVariable)
include(AddPackageDependency)

if (CLM_LIBRARIES AND CLM_INCLUDE_DIRS)

  # Do nothing. Variables are set. No need to search again

elseif (CLM_DIR)

  set(CLM_INCLUDE_DIR ${CLM_DIR}/include)
  set(CLM_LIBRARY_DIR ${CLM_DIR}/lib)
  set(CLM_TARGET clm)

  find_library(_CLM_LIBRARY
               NAMES clm
               PATHS ${CLM_DIR}/lib)

  if (_CLM_LIBRARY)
    add_imported_library(${CLM_TARGET}
                         LOCATION ${_CLM_LIBRARY}
                         LINK_LANGUAGES "Fortran")
    set(CLM_LIBRARY ${CLM_TARGET})
  endif()    
    
  # Define the LIBRARIES and INCLUDE_DIRS
  set(CLM_INCLUDE_DIRS ${CLM_INCLUDE_DIR})
  set(CLM_LIBRARIES ${CLM_LIBRARY})

endif (CLM_LIBRARIES AND CLM_INCLUDE_DIRS)

# Send useful message if everything is found
find_package_handle_standard_args(CLM DEFAULT_MSG
                                  CLM_INCLUDE_DIRS
                                  CLM_LIBRARIES)

# find_package should set CLM_FOUND but it does not!
if (CLM_LIBRARIES AND CLM_INCLUDE_DIRS)
  set(CLM_FOUND TRUE)
else()
  set(CLM_FOUND FALSE)
endif()

mark_as_advanced(
  CLM_INCLUDE_DIR
  CLM_LIBRARY_DIR
  CLM_LIBRARIES
)
