# -*- mode: cmake -*-

#
#  Amanzi
#    Chemistry
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)
include(InstallManager)

#
# Define a project name
# After this command the following varaibles are defined
#   CHEM_SOURCE_DIR
#   CHEM_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(CHEM)

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${STATE_SOURCE_DIR})

# Chemistry includes directory
# The set will define the variable in the chemistry scope
# The global property will allow other directories (MPC) to
# access this directory name through a variable.
set(CHEM_INCLUDES_DIR "${CHEM_SOURCE_DIR}/includes")
set_property(GLOBAL PROPERTY CHEM_INCLUDES_DIR "${CHEM_INCLUDES_DIR}")
include_directories(${CHEM_INCLUDES_DIR})

set(GEOCHEM_LIBS geochemrxns geochembase geochemsolvers geochemutil)
set_property(GLOBAL PROPERTY GEOCHEM_LIBS "${GEOCHEM_LIBS}")

add_subdirectory(utility)
add_subdirectory(solvers)
add_subdirectory(base_chemistry)
add_subdirectory(reactions)

#
# Install targets
#
file(GLOB chem_inc "${CHEM_INCLUDES_DIR}/*")
add_install_include_file(${chem_inc})

if (BUILD_TESTS)
  add_subdirectory(tests_unit)
  add_subdirectory(tests_regression)
endif()
