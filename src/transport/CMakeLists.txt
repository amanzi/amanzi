# -*- mode: cmake -*-

#
#  Amanzi
#   Transport process kernel
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)

#
# Define a project name
# After this command the following varaibles are defined
#   TRANSPORT_SOURCE_DIR
#   TRANSPORT_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(TRANSPORT)

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${MESH_BASE_SOURCE_DIR})
include_directories(${MESH_DATA_SOURCE_DIR})
include_directories(${STATE_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Epetra_INCLUDE_DIRS})

# Need to move this to the root level -- lpritch
add_definitions("-DUSE_MPI")

#
# Library: transport
#
add_library(transport Transport_PK.cpp Transport_State.cpp)
target_link_libraries(transport  mesh_data error_handling state flow chemistry ${Epetra_LIBRARIES})

if ( BUILD_TESTS )

    # Add UnitTest includes
    include_directories(${UnitTest_INCLUDE_DIRS})

    if (NOT (${TRANSPORT_SOURCE_DIR} EQUAL ${TRANSPORT_BINARY_DIR}) )
        execute_process(COMMAND ${CMAKE_COMMAND} -E 
          copy_directory ${TRANSPORT_SOURCE_DIR}/test ${TRANSPORT_BINARY_DIR}/test) 
    endif()

    if ( ENABLE_MOAB_Mesh )

        # Add MOAB include directory
        include_directories(${MOAB_INCLUDE_DIRS})
        include_directories(${MOAB_SOURCE_DIR})
        include_directories(${TRANSPORT_SOURCE_DIR})

        set(amanzi_libs simple_mesh moab_mesh state flow transport)

        # Test: transport
        add_executable(test_serial  test/Main.cc test/test_transport.cpp test/test_moab.cpp)
        target_link_libraries(test_serial ${amanzi_libs} ${UnitTest_LIBRARIES})
        add_unit_test(transport test_serial)

        # Test: transport (parallel)
        add_executable(test_parallel test/Main.cc test/test_moab_parallel.cpp)
        target_link_libraries(test_parallel ${amanzi_libs} ${UnitTest_LIBRARIES})
        add_unit_test(transport_parallel test_parallel NPROCS ${MPI_EXEC_MAX_NPROCS})

    endif()    

endif()


