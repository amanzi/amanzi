# -*- mode: cmake -*-

#
#  Amanzi
#   Energy process kernel
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)

#
# Define a project name
# After this command the following varaibles are defined
#   ENERGY_SOURCE_DIR
#   ENERGY_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(ENERGY)

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${MESH_DATA_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${DATA_STRUCTURES_SOURCE_DIR})
include_directories(${STATE_SOURCE_DIR})
include_directories(${WHETSTONE_SOURCE_DIR})
include_directories(${SOLVERS_SOURCE_DIR})
include_directories(${TIME_INTEGRATION_SOURCE_DIR})
include_directories(${FUNCS_SOURCE_DIR})
include_directories(${MFUNCS_SOURCE_DIR})
include_directories(${OUTPUT_SOURCE_DIR})
include_directories(${DBG_SOURCE_DIR})
include_directories(${OPERATORS_SOURCE_DIR})
include_directories(${PKS_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${ASCEMIO_INCLUDE_DIR})

#
# Library: energy
#
set(energy_src_files Energy_PK.cc)
set(energy_tpl_libs ${Teuchos_LIBRARIES} ${Epetra_LIBRARIES})
if ( Ifpack_ENABLE_HYPRE )
  list(APPEND energy_tpl_libs ${Ifpack_LIBRARIES})
endif()
add_amanzi_library(energy SOURCE ${energy_src_files}
                   LINK_LIBS functions mesh_functions ${energy_tpl_libs} time_integration data_structures state)

#
# Install Targets
#
set(flow_inc_files 
  Energy_PK.hh)

add_install_include_file(${energy_inc_files})

if (BUILD_TESTS) 
  
#if (0)
    # Add UnitTest include directoy
    include_directories(${UnitTest_INCLUDE_DIRS})
    include_directories(${MESH_FACTORY_SOURCE_DIR})

    # Copy test directory files if an out of source build
    if (NOT (${ENERGY_SOURCE_DIR} EQUAL ${ENERGY_BINARY_DIR}) )
        execute_process(COMMAND ${CMAKE_COMMAND} -E 
          copy_directory ${ENERGY_SOURCE_DIR}/test ${ENERGY_BINARY_DIR}/test) 
    endif()

    # Add the flow directory to the include paths
    include_directories(${ENERGY_SOURCE_DIR})

    set(amanzi_libs geometry mesh mesh_audit mesh_factory state energy operators whetstone data_structures)

    # Test: 1
    #add_amanzi_test(energy_test energy_test
    #                KIND int
    #                SOURCE test/Main.cc test/energy_test.cc
    #                LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})
endif()

