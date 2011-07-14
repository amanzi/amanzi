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
include_directories(${GEOMETRY_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Epetra_INCLUDE_DIRS})

# Need to move this to the root level -- lpritch
add_definitions("-DUSE_MPI" "-DMSTK_HAVE_MPI")

#
# Library: transport
#
add_library(transport Transport_PK.cpp Transport_State.cpp Reconstruction.cpp)
target_link_libraries(transport geometry mesh mesh_data error_handling state ${Epetra_LIBRARIES})

#
# Install Targets
#
file(GLOB transport_inc_files "*.hpp")
add_install_include_file(${transport_inc_files})
add_install_library(transport)

#
# Install Targets
#
file(GLOB transport_inc_files "*.hpp")
add_install_include_file(${transport_inc_files})
add_install_library(transport)

if ( BUILD_TESTS )
    # Add UnitTest includes
    include_directories(${UnitTest_INCLUDE_DIRS})

    if (NOT (${TRANSPORT_SOURCE_DIR} EQUAL ${TRANSPORT_BINARY_DIR}) )
        execute_process(COMMAND ${CMAKE_COMMAND} -E 
          copy_directory ${TRANSPORT_SOURCE_DIR}/test ${TRANSPORT_BINARY_DIR}/test) 
    endif()

    include_directories(${TRANSPORT_SOURCE_DIR})
    include_directories(${Boost_INCLUDE_DIRS})
    include_directories(${MESH_SIMPLE_SOURCE_DIR})

    set(amanzi_libs geometry mesh mesh_audit simple_mesh state transport)

    # Test 1: transport with a simple mesh
    add_executable(transport_serial test/Main.cc test/test_transport.cpp)
    target_link_libraries(transport_serial ${amanzi_libs} ${UnitTest_LIBRARIES})
    add_amanzi_test(transport_serial transport_serial KIND unit)

    # Test 2: high-order transport with a simple mesh
    add_executable(transport_high_order test/Main.cc test/test_high_order.cpp)
    target_link_libraries(transport_high_order ${amanzi_libs} ${UnitTest_LIBRARIES})
    add_amanzi_test(transport_high_order transport_high_order KIND unit)

    # Tests 3 & 4
    if ( ENABLE_MSTK_Mesh )
        include_directories(${MSTK_INCLUDE_DIRS})
        include_directories(${MSTK_SOURCE_DIR})

        set(amanzi_libs geometry mesh mesh_audit mstk_mesh state transport)

        # Test: transport (serial)
        add_executable(transport_serial_mstk test/Main.cc test/test_mstk.cpp)
        target_link_libraries(transport_serial_mstk ${amanzi_libs} ${UnitTest_LIBRARIES})
        add_amanzi_test(transport_serial_mstk transport_serial_mstk KIND unit)

        # Test: transport (parallel)
        add_executable(transport_parallel_mstk test/Main.cc test/test_mstk_parallel.cpp)
        target_link_libraries(transport_parallel_mstk ${amanzi_libs} ${UnitTest_LIBRARIES})
        add_amanzi_test(transport_parallel_mstk transport_parallel_mstk NPROCS ${MPI_EXEC_MAX_NPROCS} KIND int)
    endif()    

    # Tests 5 and 6
    if ( ENABLE_MOAB_Mesh )
        include_directories(${MOAB_INCLUDE_DIRS})
        include_directories(${MOAB_SOURCE_DIR})

        set(amanzi_libs geometry mesh mesh_audit moab_mesh state transport)

        # Test: transport (serial)
        add_executable(transport_serial_moab test/Main.cc test/test_moab.cpp)
        target_link_libraries(transport_serial_moab ${amanzi_libs} ${UnitTest_LIBRARIES})
        add_amanzi_test(transport_serial_moab transport_serial_moab KIND unit)

        # Test: transport (parallel)
        add_executable(transport_parallel_moab test/Main.cc test/test_moab_parallel.cpp)
        target_link_libraries(transport_parallel_moab ${amanzi_libs} ${UnitTest_LIBRARIES})
        add_amanzi_test(transport_parallel_moab transport_parallel_moab NPROCS ${MPI_EXEC_MAX_NPROCS} KIND int)
    endif()    

endif()


