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
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${HDF5_INCLUDE_DIRS})

#
# Library: pks
#
set(pks_src_files PK_Factory.cc
                  PK_Utils.cc
                  PK_Physical.cc)

add_amanzi_library(pks
                   SOURCE ${pks_src_files}
                   HEADERS ${pks_inc_files}
		   LINK_LIBS state data_structures error_handling atk
                             ${Epetra_LIBRARIES})

#
# Install Targets
#
set(pks_inc_files 
  CommonDefs.hh
  PK.hh
  PK_BDF.hh
  PK_Implicit.hh
  PK_Explicit.hh
  PK_Physical.hh
  PK_PhysicalExplicit.hh
  PK_PhysicalBDF.hh
  PK_DomainFunction.hh
  PK_DomainFunctionFactory.hh
  PK_DomainFunctionSimple.hh
  PK_DomainFunctionVolume.hh
  PK_DomainFunctionVolumeFraction.hh
  PK_DomainFunctionWeight.hh
  PK_DomainFunctionCoupling.hh
  PK_DomainFunctionSubgrid.hh
  PK_DomainFunctionSubgridReturn.hh
  PK_DomainFunctionFirstOrderExchange.hh
  PK_DomainFunctionSimpleWell.hh
  PK_Factory.hh
  PK_Utils.hh)

add_install_include_file(${pks_inc_files})

#
# No tests so far
#
if (BUILD_TESTS) 
endif()

if (ENABLE_AmanziPhysicsModule)
  add_subdirectory(eos)
  add_subdirectory(flow)
  add_subdirectory(multiphase)
  add_subdirectory(navier_stokes)
endif()

# chemistry isn't REALLY physics...
add_subdirectory(chemistry)

# ok, transport is physics but oh well...
add_subdirectory(transport)

# add shallow water
add_subdirectory(shallow_water)

if (ENABLE_AmanziPhysicsModule)
  add_subdirectory(energy)
  add_subdirectory(dummy)
endif()

add_subdirectory(mpc_pk)

