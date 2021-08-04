# -*- mode: cmake -*-

#
#  Amanzi
#    MOAB Mesh frameworks
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)

#
# Define a project name
# After this command the following varaibles are defined
#   MOAB_SOURCE_DIR
#   MOAB_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(MESH_MOAB)

# Amanzi include directories
include_directories(${ATK_SOURCE_DIR})
include_directories(${DBC_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${MOAB_INCLUDE_DIRS})
include_directories(${Teuchos_INCLUDE_DIRS})

# Need this define. Errors from MOAB include files 
# about MPI_COOM_WORLD. --lpritch
add_definitions("-DUSE_MPI")

#
# Library: mesh_moab
#
# The MOAB_LIBRARIES variable does not have HDF5
# libraries listed. Ticket #335
file(GLOB mesh_moab_src_files "*.cc") 
file(GLOB mesh_moab_inc_files "*.hh")
add_amanzi_library(mesh_moab
                   SOURCE ${mesh_moab_src_files} HEADERS ${mesh_moab_inc_files}
                   LINK_LIBS geometry mesh atk
                             ${MOAB_LIBRARIES} ${Epetra_LIBRARIES} ${Teuchos_LIBRARIES}) 

