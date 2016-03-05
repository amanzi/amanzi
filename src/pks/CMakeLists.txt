# -*- mode: cmake -*-

#
#  Amanzi
#    PK class
#

#
# Define a project name
# After this command the following varaibles are defined
#   PKS_SOURCE_DIR
#   PKS_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(PKS)

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${MESH_DATA_SOURCE_DIR})
include_directories(${DATA_STRUCTURES_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${WHETSTONE_SOURCE_DIR})
include_directories(${OPERATORS_SOURCE_DIR})
include_directories(${TIME_INTEGRATION_SOURCE_DIR})
include_directories(${SOLVERS_SOURCE_DIR})
include_directories(${FUNCS_SOURCE_DIR})
include_directories(${MFUNCS_SOURCE_DIR})
include_directories(${STATE_SOURCE_DIR})
include_directories(${OUTPUT_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${ASCEMIO_INCLUDE_DIR})

#
# Library: pks
#
set(pks_src_files PK_Factory.cc
                  PK_BoundaryFunction.cc
                  PK_DomainFunction.cc
                  PK_Utils.cc
                  PK_Default.cc
#                  PK_BDF.cc
                  PK_Physical.cc
#                  PK_PhysicalBDF.cc
                  PK_PhysicalBDF_ATS.cc
   )

add_amanzi_library(pks
                   SOURCE ${pks_src_files} HEADERS ${pks_inc_files}
                   LINK_LIBS )

#
# Install Targets
#
set(pks_inc_files 
  CommonDefs.hh
  PK.hh
  PK_Default.hh
  PK_BDF.hh
  PK_Implicit.hh
  PK_Explicit.hh
  PK_Physical.hh
  PK_PhysicalBDF.hh
  PK_PhysicalBDF_ATS.hh
  PK_BoundaryFunction.hh
  PK_DomainFunction.hh
  PK_Factory.hh
  PK_Utils.hh)

add_install_include_file(${pks_inc_files})

if (BUILD_TESTS) 
endif()

add_subdirectory(eos)

add_subdirectory(flow)
add_subdirectory(chemistry)
add_subdirectory(transport)
add_subdirectory(energy)
add_subdirectory(dummy)
add_subdirectory(mpc_pk)







