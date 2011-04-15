# -*- mode: cmake -*-

#
#  Amanzi
#    MSTK Mesh frameworks
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)

#
# Define a project name
# After this command the following varaibles are defined
#   MSTK_SOURCE_DIR
#   MSTK_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(MSTK)

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${MESH_BASE_SOURCE_DIR})
include_directories(${MESH_DATA_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${MSTK_INCLUDE_DIRS})

# Need this define. Errors from MSTK include files 
# about MPI_COMM_WORLD. --lpritch
add_definitions("-DUSE_MPI -DMSTK_HAVE_MPI")

#
# Library: mstk_mesh
#
set(mstk_mesh_files Cell_topology.cc Element_category.cc Mesh_maps_mstk.cc) 
add_library(mstk_mesh ${mstk_mesh_files})
target_link_libraries(mstk_mesh mesh_data error_handling ${MSTK_LIBRARIES} ${ExodusII_LIBRARIES} ${Epetra_LIBRARIES} ${Teuchos_LIBRARIES})

if (BUILD_TESTS)
    
    # Add UnitTest ExodusII include directories
    include_directories(${UnitTest_INCLUDE_DIRS})

    # Need to copy files for the tests 
    if (NOT (${MSTK_SOURCE_DIR} EQUAL ${MSTK_BINARY_DIR}))
        execute_process(COMMAND ${CMAKE_COMMAND} -E 
          copy_directory ${MSTK_SOURCE_DIR}/test ${MSTK_BINARY_DIR}/test) 
    endif()


    # Test: mstk_mesh
    add_executable(test_mstk_mesh 
                   test/Main.cc 
                   test/test_hex_2x2x2.cc
                   test/test_hex_3x3x2.cc
                   test/test_hex_4x4x4.cc)
    target_link_libraries(test_mstk_mesh mstk_mesh ${UnitTest_LIBRARIES})
    add_unit_test(mstk_mesh test_mstk_mesh )

    # Test: mstk_mesh_parallel
    add_executable(test_parallel_mstk_mesh
        test/Main.cc
        test/test_hex_4x4x4_4P.cc)
    target_link_libraries(test_parallel_mstk_mesh mstk_mesh ${UnitTest_LIBRARIES})
    add_integration_test(mstk_mesh_parallel test_parallel_mstk_mesh NPROCS 4)

endif()

