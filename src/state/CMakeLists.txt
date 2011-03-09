# -*- mode: cmake -*-

#
#  Amanzi
#    State class
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(AddParallelTest)

#
# Define a project name
# After this command the following varaibles are defined
#   STATE_SOURCE_DIR
#   STATE_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(STATE)

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${OUTPUT_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})

# Need to move this up to the root -- lpritch
add_definitions("-DUSE_MPI")


# This state library will move to a new location 
# once Markus has created the data manager.
# I use global properties since all the PROJECT_NAME
# variables in this directory are tied to the MPC.
# These properties will be replaced with the appropriate
# directory project name once it is moved.
#
# Library: state
#
add_library(state State.cpp)
set(state_tpl_list ${Boost_LIBRARIES} ${Teuchos_LIBRARIES} ${Epetra_LIBRARIES})
target_link_libraries(state output simple_mesh mesh_data cell_geometry ${state_tpl_libs})
set_property(GLOBAL PROPERTY STATE_SOURCE_DIR ${MPC_SOURCE_DIR})
set_property(GLOBAL PROPERTY STATE_BINARY_DIR ${MPC_BINARY_DIR})
