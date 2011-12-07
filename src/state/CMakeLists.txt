# -*- mode: cmake -*-

#
#  ATS
#    State class
#

#
# Define a project name
# After this command the following varaibles are defined
#   STATE_SOURCE_DIR
#   STATE_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(STATE)

# ATS include directories
# none required

# Amanzi install include path
include_directories(${Amanzi_INCLUDE_DIR})

# External (TPL) include directories
include_directories(${Amanzi_TPL_Teuchos_INCLUDE_DIRS})
include_directories(${Amanzi_TPL_Epetra_INCLUDE_DIRS})
include_directories(${Amanzi_TPL_HDF5_INCLUDE_DIRS})
include_directories(${Amanzi_TPL_ASCEMIO_INCLUDE_DIR})

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
add_library(state Field.cc State.cc Vis.cc)
add_install_library(state)
set(state_tpl_list ${Amanzi_TPL_Boost_LIBRARIES} ${Amanzi_TPL_Teuchos_LIBRARIES} ${Amanzi_TPL_Epetra_LIBRARIES})
target_link_libraries(state output simple_mesh mesh_data cell_geometry ${state_tpl_libs})
set_property(GLOBAL PROPERTY STATE_SOURCE_DIR ${MPC_SOURCE_DIR})
set_property(GLOBAL PROPERTY STATE_BINARY_DIR ${MPC_BINARY_DIR})

#
# Install Target
#
file(GLOB state_inc_files "*.hh")
add_install_include_file(${state_inc_files})

