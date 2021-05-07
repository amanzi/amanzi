# -*- mode: cmake -*-

#
# Amanzi CrunchTope Find Module
#
# Usage:
#    Control the search through CrunchTope_DIR installation prefix.
#
#    This module does not search default paths! 
#
#    Following variables are set:
#    CrunchTope_FOUND          (BOOL)    Flag indicating if CrunchTope was found
#    CrunchTope_INCLUDE_DIR    (PATH)    Path to the CrunchTope include file
#    CrunchTope_INCLUDE_DIRS   (LIST)    List of all required include files
#    CrunchTope_LIBRARY_DIR    (PATH)    Path to the CrunchTope library
#    CrunchTope_LIBRARIES      (LIST)    List of all required CrunchTope libraries
#
# #############################################################################

# Standard CMake modules see CMAKE_ROOT/Modules
include(FindPackageHandleStandardArgs)

# Amanzi CMake functions see <root>/tools/cmake for source
include(PrintVariable)
include(AddPackageDependency)

if (CrunchTope_LIBRARIES AND CrunchTope_INCLUDE_DIRS)

  # Do nothing. Variables are set. No need to search again

elseif (CrunchTope_DIR)

  set(CrunchTope_INCLUDE_DIR ${CrunchTope_DIR}/include)
  set(CrunchTope_LIBRARY_DIR ${CrunchTope_DIR}/lib)
  set(CrunchTope_TARGET crunchchem)

  find_library(_CrunchTope_LIBRARY
               NAMES crunchchem
               PATHS ${CrunchTope_DIR}/lib)

  if (_CrunchTope_LIBRARY)
    add_imported_library(${CrunchTope_TARGET}
                         LOCATION ${_CrunchTope_LIBRARY}
                         LINK_LANGUAGES "Fortran")
    set(CrunchTope_LIBRARY ${CrunchTope_TARGET})
  endif()    
    
  # Define the LIBRARIES and INCLUDE_DIRS
  set(CrunchTope_INCLUDE_DIRS ${CrunchTope_INCLUDE_DIR})
  set(CrunchTope_LIBRARIES ${CrunchTope_LIBRARY})

endif (CrunchTope_LIBRARIES AND CrunchTope_INCLUDE_DIRS)

# Send useful message if everything is found
find_package_handle_standard_args(CrunchTope DEFAULT_MSG
                                  CrunchTope_INCLUDE_DIRS
                                  CrunchTope_LIBRARIES)

# find_package should set CrunchTope_FOUND but it does not!
if (CrunchTope_LIBRARIES AND CrunchTope_INCLUDE_DIRS)
  set(CrunchTope_FOUND TRUE)
else()
  set(CrunchTope_FOUND FALSE)
endif()

mark_as_advanced(
  CrunchTope_INCLUDE_DIR
  CrunchTope_LIBRARY_DIR
  CrunchTope_LIBRARIES
)
