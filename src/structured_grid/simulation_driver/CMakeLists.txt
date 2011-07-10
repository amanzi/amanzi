# -*- mode: cmake -*-
#
# Define a project name
# After this command the following varaibles are defined
#   BOXDRIV_SOURCE_DIR
#   BOXDRIV_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(BOXDRIV)

get_property(CHEM_INCLUDES_DIR GLOBAL PROPERTY CHEM_INCLUDES_DIR)
include_directories(${CHEM_INCLUDES_DIR})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CCSE_INCLUDE_DIRS})

include_directories(${CHEM_INCLUDES_DIR})
include_directories(${PMAMR_INCLUDE_DIRS})

include_directories(${PLATFORM_INTERFACE_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})

add_library(boxdriv amanzi_structured_grid_simulation_driver.cpp)

if (BUILD_TESTS)

endif()

