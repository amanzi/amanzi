# -*- mode: cmake -*-

#
#  Amanzi
#    PK class
#

# Define a project name
# After this command the following varaibles are defined
#   PKS_SOURCE_DIR
#   PKS_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(PKS)

# Amanzi include directories
include_directories(${ATK_SOURCE_DIR})
include_directories(${DATA_STRUCTURES_SOURCE_DIR})
include_directories(${DBC_SOURCE_DIR})
include_directories(${FUNCS_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${MESH_DATA_SOURCE_DIR})
include_directories(${MFUNCS_SOURCE_DIR})
include_directories(${OPERATORS_SOURCE_DIR})
include_directories(${OUTPUT_SOURCE_DIR})
include_directories(${SOLVERS_SOURCE_DIR})
include_directories(${STATE_SOURCE_DIR})
include_directories(${STATE_SOURCE_DIR}/data)
include_directories(${STATE_SOURCE_DIR}/evaluators)
include_directories(${TIME_INTEGRATION_SOURCE_DIR})
include_directories(${WHETSTONE_SOURCE_DIR})

# External (TPL) include directories
include_directories(${ASCEMIO_INCLUDE_DIR})
include_directories(${Tpetra_INCLUDE_DIRS})
include_directories(${HDF5_INCLUDE_DIRS})

#
# Library: pks
#
set(pks_src_files PKFactory.cc)

add_amanzi_library(pks
                   SOURCE ${pks_src_files}
                   HEADERS ${pks_inc_files}
		   LINK_LIBS state data_structures error_handling atk
                             ${Tpetra_LIBRARIES})

#
# Install Targets
#
set(pks_inc_files 
  PK.hh
  PK_BDF.hh
  PK_Implicit.hh
  PK_Explicit.hh
  PK_Physical.hh
  PK_PhysicalExplicit.hh
  PK_PhysicalBDF.hh
  PKFactory.hh
)

add_install_include_file(${pks_inc_files})


