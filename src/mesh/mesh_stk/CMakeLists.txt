# -*- mode: cmake -*-
#
#  Amanzi
#    STK Mesh framework (Trilinos package)
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(TestManager)
include(PrintVariable)

#
# Define a project name
# After this command the following varaibles are defined
#   STK_SOURCE_DIR
#   STK_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(STK)

# Amanzi include directories
include_directories(${ATK_SOURCE_DIR})
include_directories(${DBC_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${MESH_BASE_SOURCE_DIR})
include_directories(${MESH_DATA_SOURCE_DIR})
include_directories(${EXO_SOURCE_DIR})


# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${STK_INCLUDE_DIRS})
include_directories(${ExodusII_INCLUDE_DIRS})

#
# Library: stk_mesh
#
file(GLOB stk_mesh_files "*.cc") 
file(GLOB stk_inc_files "*.hh")
add_amanzi_library(stk_mesh
                   SOURCE ${stk_mesh_files} HEADERS ${stk_inc_files}
                   LINK_LIBS geometry mesh mesh_data ${Teuchos_LIBRARIES} ${Epetra_LIBRARIES} ${STK_LIBRARIES} ${ExodusII_LIBRARIES})

if (BUILD_TESTS)

    # Add UnitTest includes
    include_directories(${UnitTest_INCLUDE_DIRS})

   # Need to copy files for the tests 
    if (NOT (${STK_SOURCE_DIR} EQUAL ${STK_BINARY_DIR}))
        execute_process(COMMAND ${CMAKE_COMMAND} -E 
          copy_directory ${STK_SOURCE_DIR}/test ${STK_BINARY_DIR}/test) 
    endif()



    # Test: stk_mesh
    add_amanzi_test(stk_mesh test_stk_mesh
                    KIND unit
                    SOURCE
                          test/Main.cc
                          test/test_Mesh.cc
                          test/test_Maps.cc
                          test/test_Hex.cc
                          test/test_STK.cc
                    LINK_LIBS
                          stk_mesh mesh_audit
                          ${UnitTest_LIBRARIES} )

#  # parallel tests are limited to 4 processors because the read tests
#  # use files with a maximum of 4 partitions and the internally
#  # generated meshes are small (8 cells)
#
#  set(LOCAL_MPI_EXEC_MAX_NUMPROCS MPI_EXEC_MAX_NUMPROCS)
#  if (LOCAL_MPI_EXEC_MAX_NUMPROCS GREATER 4)
#    set(LOCAL_MPI_EXEC_MAX_NUMPROCS 4)
#  endif()

    # Test: stk_mesh_parallel
    add_amanzi_test(stk_mesh_parallel test_stk_mesh_parallel
                    KIND unit
                    NPROCS 4
                    SOURCE test/Main.cc test/test_Hex.cc
                    LINK_LIBS stk_mesh mesh_audit
                              ${UnitTest_LIBRARIES})
     
    # Test: stk_mesh_read (serial and parallel)
    add_amanzi_test(stk_mesh_read test_stk_mesh_read
                    KIND int
                    SOURCE test/Main.cc test/test_Read.cc
                    LINK_LIBS exodus_reader stk_mesh mesh_audit ${UnitTest_LIBRARIES})
    add_amanzi_test(stk_mesh_read_parallel test_stk_mesh_read NPROCS 4 KIND int)

    # Test: Create sets and respond to set queries
    add_amanzi_test(stk_mesh_sets test_stk_sets
                    KIND int
                    SOURCE test/Main.cc test/test_Sets.cc
                    LINK_LIBS stk_mesh exodus_reader ${UnitTest_LIBRARIES})

endif()

