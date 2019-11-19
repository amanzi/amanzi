# -*- mode: cmake -*-

#
#  
#   Transport process kernel
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
##include(PrintVariable)
##include(TestManager)

#
# Define a project name
# After this command the following varaibles are defined
#   TRANSPORT_SOURCE_DIR
#   TRANSPORT_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
##project(TRANSPORT)


# External (TPL) include directories
##include_directories(${Epetra_INCLUDE_DIRS})
##include_directories(${ASCEMIO_INCLUDE_DIR})
# Amanzi include directories
##include_directories(${ATK_SOURCE_DIR})
##include_directories(${CHEM_INCLUDES_DIR})
##include_directories(${CHEMPK_SOURCE_DIR})
##include_directories(${DATA_STRUCTURES_SOURCE_DIR})
##include_directories(${FUNCS_SOURCE_DIR})
##include_directories(${GEOMETRY_SOURCE_DIR})
##include_directories(${MESH_SOURCE_DIR})
##include_directories(${MESH_FACTORY_SOURCE_DIR})
##include_directories(${MFUNCS_SOURCE_DIR})
##include_directories(${OUTPUT_SOURCE_DIR})
##include_directories(${OPERATORS_SOURCE_DIR})
##include_directories(${PKS_SOURCE_DIR})
##include_directories(${SOLVERS_SOURCE_DIR})
##include_directories(${STATE_SOURCE_DIR})
##include_directories(${TIME_INTEGRATION_SOURCE_DIR})
##include_directories(${WHETSTONE_SOURCE_DIR})

#include_directories(${Amanzi_TPL_MSTK_INCLUDE_DIRS})

#
# Transport registrations
#



#
# Library: transport
#
set(transport_inc_files
  sediment_transport_pk.hh
  erosion_evaluator.hh settlement_evaluator.hh trapping_evaluator.hh organic_matter_evaluator.hh
  )

set(transport_src_files sediment_transport_pk.cc
  Transport_Diffusion.cc
  erosion_evaluator.cc settlement_evaluator.cc trapping_evaluator.cc organic_matter_evaluator.cc)

add_library(pk_sediment_transport 
            ${transport_src_files}
)

install(TARGETS pk_sediment_transport DESTINATION lib)

#                   LINK_LIBS geometry mesh error_handling state data_structures 
#                             ${geochem_lib} chemistry_pk time_integration ${Epetra_LIBRARIES})

if ( 0 )
#if ( BUILD_TESTS )
    # Add UnitTest includes
    include_directories(${UnitTest_INCLUDE_DIRS})

endif()
