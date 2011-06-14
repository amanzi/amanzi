# -*- mode: cmake -*-

#
#  Amanzi
#    Mesh 
#
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)

#
# Define a project name
# After this command the following varaibles are defined
#   MESH_SOURCE_DIR
#   MESH_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(MESH)

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})

# External (TPL) directories
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})

#
# Library: cell_geometry
#
add_library(cell_geometry cell_geometry.cc)

#
# Library: mesh_audit
#
add_library(mesh_audit MeshAudit.cc)
target_link_libraries(mesh_audit error_handling cell_geometry ${Epetra_LIBRARIES})

#
# Install Targets
#
file(GLOB mesh_inc_files "*.hh")
add_install_include_file(${mesh_inc_files})
add_install_library(cell_geometry mesh_audit)
