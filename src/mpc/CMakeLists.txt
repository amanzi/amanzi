# -*- mode: cmake -*-

#
#  Amanzi
#    Mulitprocess Coordinator (MPC)
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)

#
# Define a project name
# After this command the following varaibles are defined
#   MPC_SOURCE_DIR
#   MPC_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(MPC)

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${SIMPLE_MESH_SOURCE_DIR})
include_directories(${MESH_DATA_SOURCE_DIR})
include_directories(${OUTPUT_SOURCE_DIR})
include_directories(${STATE_SOURCE_DIR})
include_directories(${CHEM_SOURCE_DIR})
include_directories(${TRANSPORT_SOURCE_DIR})
include_directories(${FLOW_SOURCE_DIR})

# Chemistry includes directory
get_property(CHEM_INCLUDES_DIR GLOBAL PROPERTY CHEM_INCLUDES_DIR)
include_directories(${CHEM_INCLUDES_DIR})

# External (TPL) include directories
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
if(ENABLE_CGNS)
    include_directories(${CGNS_INCLUDES})
endif()    

# Need to move this up to the root -- lpritch
add_definitions("-DUSE_MPI")


#
# Library: mpc
#
add_library(mpc MPC.cpp)
set(mpc_tpl_libs ${Boost_LIBRARIES} ${Teuchos_LIBRARIES} ${Epetra_LIBRARIES})
if (ENABLE_CGNS)
    list(APPEND mpc_tpl_libs ${CGNS_LIBRARIES})
endif()    
list(REMOVE_DUPLICATES mpc_tpl_libs) # Concerned about command line length! --lpritch
target_link_libraries(mpc output state simple_mesh mesh_data cell_geometry ${mpc_tpl_libs})

if (BUILD_TESTS)
    
    # Add UnitTest includes
    include_directories(${UnitTest_INCLUDE_DIRS})

    # Add MPC include directory. Need to remove this required path. -- lpritch
    include_directories(${MPC_SOURCE_DIR})

    # Copy test subdirectory for out of source builds
    if (NOT (MPC_SOURCE_DIR STREQUAL MPC_BINARY_DIR))
        execute_process(COMMAND ${CMAKE_COMMAND} -E 
          copy_directory ${MPC_SOURCE_DIR}/test ${MPC_BINARY_DIR}/test) 
    endif()

    
    if ( ENABLE_MOAB_Mesh )
        include_directories(${MOAB_SOURCE_DIR})
        include_directories(${MOAB_INCLUDE_DIRS})

        # Test: mpc-moab
        add_executable(test_mpc_moab test/Main.cc test/test_driver.cpp)
        target_link_libraries(test_mpc_moab mpc moab_mesh chemistry transport flow ${UnitTest_LIBRARIES})
        add_integration_test(mpc-moab test_mpc_moab)
    
    endif()

endif()
