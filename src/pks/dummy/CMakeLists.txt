# -*- mode: cmake -*-

#
#  Amanzi
#   Dummy process kernel
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)
#include(LibraryManager)

#
# Define a project name
# After this command the following varaibles are defined
#   DUMMYPK_SOURCE_DIR
#   DUMMYPK_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(DUMMYPK)

# Amanzi include directories
include_directories(${ATK_SOURCE_DIR})
include_directories(${DBC_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${MESH_DATA_SOURCE_DIR})
include_directories(${MPC_SOURCE_DIR})
include_directories(${PKS_SOURCE_DIR})
include_directories(${SOLVERS_SOURCE_DIR})
include_directories(${STATE_SOURCE_DIR})
include_directories(${TIME_INTEGRATION_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${HDF5_C_INCLUDE_DIR})
include_directories(${ASCEMIO_INCLUDE_DIR})

#
# Dummy registrations
#
register_evaluator_with_factory(HEADERFILE Dummy_PK_reg.hh
  LISTNAME   PKS_DUMMY_REGISTRATIONS)

generate_evaluators_registration_header( 
  HEADERFILE pks_dummy_registration.hh   
  LISTNAME   PKS_DUMMY_REGISTRATIONS
  INSTALL    True )


#
# Library: dummy
#
set(dummy_src_files Dummy_PK.cc)
add_amanzi_library(dummy 
                   SOURCE ${dummy_src_files}
		   LINK_LIBS pks operators ${Epetra_LIBRARIES})
#
# Install Targets
#
set(dummy_inc_files Dummy_PK.hh)  

add_install_include_file(${dummy_inc_files})


