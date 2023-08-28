# -*- mode: cmake -*-
#
# Define a project name
# After this command the following varaibles are defined
#   SDRIVER_SOURCE_DIR
#   SDRIVER_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(SDRIVER)

get_property(CHEM_INCLUDES_DIR GLOBAL PROPERTY CHEM_INCLUDES_DIR)
include_directories(${CHEM_INCLUDES_DIR})
include_directories(${CHEMPK_SOURCE_DIR})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CCSE_INCLUDE_DIRS})

include_directories(${CHEM_INCLUDES_DIR})
include_directories(${PMAMR_INCLUDE_DIRS})
include_directories(${PLATFORM_INTERFACE_SOURCE_DIR}) # havent yet visited these folders
include_directories(${ATK_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Trilinos_INCLUDE_DIRS})
include_directories(${PETSc_INCLUDE_DIRS})

set(CXX_header_files amanzi_structured_grid_simulation_driver.H)
set(CXX_source_files amanzi_structured_grid_simulation_driver.cpp)

add_library(sdriver OBJECT ${CXX_source_files}
                           ${CXX_headers_files})

if (BUILD_TESTS)

endif()

