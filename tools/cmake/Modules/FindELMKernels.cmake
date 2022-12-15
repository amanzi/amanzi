# -*- mode: cmake -*-

#
# Amanzi ELMKernels Find Module
#
# Usage:
#    Control the search through ELMKERNELS_DIR installation prefix.
#
#    This module does not search default paths! 
#
#    Following variables are set:
#    ELMKERNELS_FOUND          (BOOL)    Flag indicating if ELM was found
#    ELMKERNELS_INCLUDE_DIR    (PATH)    Path to the ELM include file
#    ELMKERNELS_INCLUDE_DIRS   (LIST)    List of all required include files
#    ELMKERNELS_LIBRARY_DIR    (PATH)    Path to the ELM library
#    ELMKERNELS_LIBRARIES      (LIST)    List of all required ELM libraries
#
# #############################################################################

# Standard CMake modules see CMAKE_ROOT/Modules
include(FindPackageHandleStandardArgs)

# Amanzi CMake functions see <root>/tools/cmake for source
include(PrintVariable)
include(AddPackageDependency)

if (ELMKERNELS_LIBRARIES AND ELMKERNELS_INCLUDE_DIRS)

  # Do nothing. Variables are set. No need to search again

elseif (ELMKERNELS_DIR)

  set(ELMKERNELS_INCLUDE_DIR ${ELMKERNELS_DIR}/include)
  set(ELMKERNELS_LIBRARY_DIR ${ELMKERNELS_DIR}/lib)
  set(ELMKERNELS_TARGET ELMKernels)

  find_library(_ELMKERNELS_LIBRARY
               NAMES ELMKernels
               PATHS ${ELMKERNELS_DIR}/lib)

  if (_ELMKERNELS_LIBRARY)
    add_imported_library(${ELMKERNELS_TARGET}
                         LOCATION ${_ELMKERNELS_LIBRARY}
                         LINK_LANGUAGES "Fortran")
    set(ELMKERNELS_LIBRARY ${ELMKERNELS_TARGET})
  endif()    
    
  # Define the LIBRARIES and INCLUDE_DIRS
  set(ELMKERNELS_INCLUDE_DIRS ${ELMKERNELS_INCLUDE_DIR})
  set(ELMKERNELS_LIBRARIES ${ELMKERNELS_LIBRARY})

endif (ELMKERNELS_LIBRARIES AND ELMKERNELS_INCLUDE_DIRS)

# Send useful message if everything is found
find_package_handle_standard_args(ELMKERNELS DEFAULT_MSG
                                  ELMKERNELS_INCLUDE_DIRS
                                  ELMKERNELS_LIBRARIES)

# find_package should set ELMKERNELS_FOUND but it does not!
if (ELMKERNELS_LIBRARIES AND ELMKERNELS_INCLUDE_DIRS)
  set(ELMKERNELS_FOUND TRUE)
else()
  set(ELMKERNELS_FOUND FALSE)
endif()

mark_as_advanced(
  ELMKERNELS_INCLUDE_DIR
  ELMKERNELS_LIBRARY_DIR
  ELMKERNELS_LIBRARIES
)
