# -*- mode: cmake -*-

#
#  Amanzi
#    Mesh Base Class
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)


#
# Define a project name
# After this command the following varaibles are defined
#   MESH_SIMPLE_SOURCE_DIR
#   MESH_SIMPLE_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(MESH_SIMPLE)

# Amanzi include directories
include_directories(${ATK_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})

# Library: simple_mesh
file(GLOB mesh_simple_source_files "*.cc")
file(GLOB mesh_simple_inc_files "*.hh")
add_amanzi_library(mesh_simple
                   SOURCE ${mesh_simple_source_files} HEADERS ${mesh_simple_inc_files}
                   LINK_LIBS mesh geometry error_handling atk 
                             ${Epetra_LIBRARIES} ${Teuchos_LIBRARIES})
