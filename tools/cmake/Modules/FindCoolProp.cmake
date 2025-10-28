# -*- mode: cmake -*-

#
# Amanzi CoolProp Find Module
#
# Usage:
#    Control the search through COOLPROP_DIR or setting environment variable
#    COOLPROP_ROOT to the COOLPROP installation prefix.
#
#    This module does not search default paths! 
#
#    Following variables are set:
#    COOLPROP_FOUND            (BOOL)       Flag indicating if COOLPROP was found
#    COOLPROP_INCLUDE_DIR      (PATH)       Path to the COOLPROP include file
#    COOLPROP_INCLUDE_DIRS     (LIST)       List of all required include files
#    COOLPROP_LIBRARY_DIR      (PATH)       Path to the COOLPROP library
#    COOLPROP_LIBRARIES        (LIST)       List of all required COOLPROP libraries
#
# #############################################################################

# Standard CMake modules see CMAKE_ROOT/Modules
include(FindPackageHandleStandardArgs)

# Amanzi CMake functions see <root>/tools/cmake for source
include(PrintVariable)
include(AddPackageDependency)

if (COOLPROP_LIBRARIES AND COOLPROP_INCLUDE_DIRS)

  # Do nothing. Variables are set. No need to search again

elseif(COOLPROP_DIR)

  set(COOLPROP_INCLUDE_DIR ${COOLPROP_DIR}/include)
  set(COOLPROP_LIBRARY_DIR ${COOLPROP_DIR}/lib)
  set(COOLPROP_TARGET CoolProp)

  find_library(_COOLPROP_LIBRARY
               NAMES CoolProp
               PATHS ${COOLPROP_LIBRARY_DIR})

  if (_COOLPROP_LIBRARY)
    set(COOLPROP_LIBRARY ${_COOLPROP_LIBRARY})
    add_imported_library(${COOLPROP_TARGET} 
                         LOCATION ${_COOLPROP_LIBRARY}
                         LINK_LANGUAGES "C;CXX")
  endif()    

  # Define the LIBRARIES and INCLUDE_DIRS
  set(COOLPROP_INCLUDE_DIRS ${COOLPROP_INCLUDE_DIR})
  set(COOLPROP_LIBRARIES ${COOLPROP_LIBRARY})

endif(COOLPROP_LIBRARIES AND COOLPROP_INCLUDE_DIRS)    

# Send useful message if everything is found
find_package_handle_standard_args(COOLPROP DEFAULT_MSG
                                  COOLPROP_INCLUDE_DIRS
                                  COOLPROP_LIBRARIES)

# find_package)handle)standard_args should set COOLPROP_FOUND but it does not!
if (COOLPROP_LIBRARIES AND COOLPROP_INCLUDE_DIRS)
  set(COOLPROP_FOUND TRUE)
else()
  set(COOLPROP_FOUND FALSE)
endif()

mark_as_advanced(
  COOLPROP_INCLUDE_DIR
  COOLPROP_LIBRARY_DIR
  COOLPROP_LIBRARY
  COOLPROP_LIBRARIES
)
