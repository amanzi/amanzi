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
# header files
#
set(pks_inc_files
  PK_Helpers.hh
  PK.hh
  PK_Default.hh
  PK_BDF.hh
  PK_BDF_Default.hh
  PK_Physical_Default.hh
  PK_Physical_Default_decl.hh
  PK_Physical_Default_impl.hh
  PK_PhysicalBDF_Default.hh
  MPC.hh
  MPCWeak.hh
  MPCStrong.hh
  PKFactory.hh
  pks_registration.hh
)


#
# implementation files
#
set(pks_src_files
  PK_Helpers.cc
  PK_BDF_Default.cc
  PK_PhysicalBDF_Default.cc
  MPCWeak.cc
  PKFactory.cc
)

add_amanzi_library(pks
                   SOURCE ${pks_src_files}
                   HEADERS ${pks_inc_files}
		   LINK_LIBS time_integration solvers state data_structures error_handling atk
                             ${Tpetra_LIBRARIES})


add_install_include_file(${pks_inc_files})


