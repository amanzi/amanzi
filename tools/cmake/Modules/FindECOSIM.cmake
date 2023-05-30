# -*- mode: cmake -*-

#
# Amanzi ECOSIM Find Module
#
# Usage:
#    Control the search through ECOSIM_DIR installation prefix.
#
#    This module does not search default paths!
#
#    Following variables are set:
#    ECOSIM_FOUND          (BOOL)    Flag indicating if ECOSIM was found
#    ECOSIM_INCLUDE_DIR    (PATH)    Path to the ECOSIM include file
#    ECOSIM_INCLUDE_DIRS   (LIST)    List of all required include files
#    ECOSIM_LIBRARY_DIR    (PATH)    Path to the ECOSIM library
#    ECOSIM_LIBRARIES      (LIST)    List of all required ECOSIM libraries
#
# #############################################################################

# Standard CMake modules see CMAKE_ROOT/Modules
include(FindPackageHandleStandardArgs)

# Amanzi CMake functions see <root>/tools/cmake for source
include(PrintVariable)
include(AddPackageDependency)

set(ECOSIM_DIR /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/local)

message("Starting FindEcoSIM...")
message("EcoSIM dir: ${ECOSIM_DIR}")

if (ECOSIM_LIBRARIES AND ECOSIM_INCLUDE_DIRS)
  message("Found EcoSIM Libs and Inc")
  # Do nothing. Variables are set. No need to search again

elseif (ECOSIM_DIR)
  message("Found EcoSIM_DIR, setting libraries and include")
  set(ECOSIM_INCLUDE_DIR ${ECOSIM_DIR}/include)
  set(ECOSIM_LIBRARY_DIR ${ECOSIM_DIR}/lib)
  set(ECOSIM_LIBRARIES ${ECOSIM_DIR}/lib)
  set(ECOSIM_TARGET ecosim)
  set(ECOSIM_INCLUDE_DIRS ${ECOSIM_DIR}/include)  

  #find_library(_ECOSIM_LIBRARY
  #             NAMES eocsim
  #            PATHS ${ECOSIM_DIR}/lib)

  message("finding library")
  find_library(_ECOSIM_LIBRARY
	       NAMES ecosim
               PATH_SUFFIXES "lib" "Lib"
               NO_DEFAULT_PATH
	       )

  message("adding imported library")	
  if (_ECOSIM_LIBRARY)
    add_imported_library(${ECOSIM_TARGET}
                         LOCATION ${_ECOSIM_LIBRARY}
                         LINK_LANGUAGES "Fortran")
    set(ECOSIM_LIBRARY ${ECOSIM_TARGET})
  endif()

  # Define the LIBRARIES and INCLUDE_DIRS
  #set(ECOSIM_INCLUDE_DIRS ${ECOSIM_INCLUDE_DIR})
  #set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARY})

endif (ECOSIM_LIBRARIES AND ECOSIM_INCLUDE_DIRS)

# Send useful message if everything is found
message("finding package handle standard args")
find_package_handle_standard_args(ECOSIM DEFAULT_MSG
                                  ECOSIM_INCLUDE_DIRS
                                  ECOSIM_LIBRARIES)

message("setting if statement")
# find_package should set ECOSIM_FOUND but it does not!
if (ECOSIM_LIBRARIES AND ECOSIM_INCLUDE_DIRS)
  set(ECOSIM_FOUND TRUE)
else()
  set(ECOSIM_FOUND FALSE)
endif()

message("running mark as advanced")

mark_as_advanced(
  ECOSIM_INCLUDE_DIR
  ECOSIM_LIBRARY_DIR
  ECOSIM_LIBRARIES
)
