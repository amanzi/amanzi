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
add_library(stk_mesh ${stk_mesh_files})
target_link_libraries(stk_mesh error_handling mesh_data mesh 
                      ${Trilinos_LIBRARIES} ${ExodusII_LIBRARIES})

#
# Install Targets
#
file(GLOB stk_inc_files "*.hh")
add_install_include_file(${stk_inc_files})
add_install_library(stk_mesh)

if (BUILD_TESTS)

    # Add UnitTest includes
    include_directories(${UnitTest_INCLUDE_DIRS})

    # Test: stk_mesh
    add_executable(test_stk_mesh test/Main.cc test/test_Mesh.cc
                     test/test_Maps.cc test/test_Hex.cc)
    target_link_libraries(test_stk_mesh stk_mesh mesh_audit
                          ${UnitTest_LIBRARIES} )
    add_amanzi_test(stk_mesh test_stk_mesh KIND unit)

#  # parallel tests are limited to 4 processors because the read tests
#  # use files with a maximum of 4 partitions and the internally
#  # generated meshes are small (8 cells)
#
#  set(LOCAL_MPI_EXEC_MAX_NUMPROCS MPI_EXEC_MAX_NUMPROCS)
#  if (LOCAL_MPI_EXEC_MAX_NUMPROCS GREATER 4)
#    set(LOCAL_MPI_EXEC_MAX_NUMPROCS 4)
#  endif()

    # Test: stk_mesh_parallel
    add_executable(test_stk_mesh_parallel test/Main.cc test/test_Hex.cc)
    target_link_libraries(test_stk_mesh_parallel stk_mesh mesh_audit
                           ${UnitTest_LIBRARIES})
    add_amanzi_test(stk_mesh_parallel test_stk_mesh_parallel NPROCS 4 KIND unit)
     
    # Test: stk_mesh_read (serial and parallel)
    add_executable(test_stk_mesh_read test/Main.cc test/test_Read.cc)
    target_link_libraries(test_stk_mesh_read exodus stk_mesh mesh_audit ${UnitTest_LIBRARIES})
    add_amanzi_test(stk_mesh_read test_stk_mesh_read KIND int)
    add_amanzi_test(stk_mesh_read_parallel test_stk_mesh_read NPROCS 4 KIND int)

endif()

