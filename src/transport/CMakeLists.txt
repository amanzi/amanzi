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
include_directories(${WHETSTONE_SOURCE_DIR})
include_directories(${EXO_SOURCE_DIR})
include_directories(${OUTPUT_SOURCE_DIR})
include_directories(${TIME_INTEGRATION_SOURCE_DIR})
include_directories(${FUNCS_SOURCE_DIR})
include_directories(${MFUNCS_SOURCE_DIR})


# External (TPL) include directories
include_directories(${Epetra_INCLUDE_DIRS})

# Need to move this to the root level -- lpritch
add_definitions("-DUSE_MPI" "-DMSTK_HAVE_MPI")

#
# Library: transport
#
add_library(transport Transport_PK.cpp Transport_State.cpp Transport_TI.cpp 
                      Transport_Limiters.cpp Transport_VandV.cpp Transport_IO.cpp
                      Reconstruction.cpp)
target_link_libraries(transport geometry mesh mesh_data error_handling state time_integration ${Epetra_LIBRARIES})
add_install_library(transport)
#
# Install Targets
#
file(GLOB transport_inc_files "*.hpp")
add_install_include_file(${transport_inc_files})

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

    set(amanzi_libs geometry mesh mesh_audit simple_mesh state transport whetstone)

    # Test: miscaleneous transport routines
    add_executable(transport_mics test/Main.cc test/transport_mics.cpp)
    target_link_libraries(transport_mics ${amanzi_libs} ${UnitTest_LIBRARIES})
    add_amanzi_test(transport_mics transport_mics KIND unit)

    add_executable(transport_advance_simple test/Main.cc test/transport_advance_simple.cpp)
    target_link_libraries(transport_advance_simple ${amanzi_libs} ${UnitTest_LIBRARIES})
    add_amanzi_test(transport_advance_simple transport_advance_simple KIND unit)

    add_executable(transport_convergence test/Main.cc test/transport_convergence.cpp)
    target_link_libraries(transport_convergence ${amanzi_libs} ${UnitTest_LIBRARIES})
    add_amanzi_test(transport_convergence transport_convergence KIND unit)

    # add_executable(transport_dispersion test/Main.cc test/transport_dispersion.cpp)
    # target_link_libraries(transport_dispersion ${amanzi_libs} ${UnitTest_LIBRARIES})
    # add_amanzi_test(transport_dispersion transport_dispersion KIND unit)

    if ( ENABLE_MSTK_Mesh )
        include_directories(${MSTK_INCLUDE_DIRS})
        include_directories(${MSTK_SOURCE_DIR})

        set(amanzi_libs geometry mesh mesh_audit mstk_mesh state transport whetstone output)

        add_executable(transport_advance_mstk test/Main.cc test/transport_advance_mstk.cpp)
        target_link_libraries(transport_advance_mstk ${amanzi_libs} ${UnitTest_LIBRARIES})
        add_amanzi_test(transport_advance_mstk transport_advance_mstk KIND unit)

        add_executable(transport_parallel_mstk test/Main.cc test/transport_parallel_mstk.cpp)
        target_link_libraries(transport_parallel_mstk ${amanzi_libs} ${UnitTest_LIBRARIES})
        add_amanzi_test(transport_parallel_mstk transport_parallel_mstk NPROCS 4 KIND int)

        add_executable(transport_2D test/Main.cc test/transport_2D.cpp)
        target_link_libraries(transport_2D ${amanzi_libs} ${UnitTest_LIBRARIES})
        add_amanzi_test(transport_2D transport_2D KIND unit)

        add_executable(transport_2D_long test/Main.cc test/transport_2D_long.cpp)
        target_link_libraries(transport_2D_long ${amanzi_libs} ${UnitTest_LIBRARIES})
        add_amanzi_test(transport_2D_long transport_2D_long KIND unit)

        add_executable(transport_3D_long test/Main.cc test/transport_3D_long.cpp)
        target_link_libraries(transport_3D_long ${amanzi_libs} ${UnitTest_LIBRARIES})
        add_amanzi_test(transport_3D_long transport_3D_long KIND unit)
    endif()    

    if ( ENABLE_STK_Mesh )
        include_directories(${STK_INCLUDE_DIRS})
        include_directories(${STK_SOURCE_DIR})

        set(amanzi_libs geometry mesh mesh_audit stk_mesh exodus_reader state transport whetstone)

        add_executable(transport_advance_stk test/Main.cc test/transport_advance_stk.cpp)
        target_link_libraries(transport_advance_stk ${amanzi_libs} ${UnitTest_LIBRARIES})
        add_amanzi_test(transport_advance_stk transport_advance_stk KIND unit)

        add_executable(transport_parallel_stk test/Main.cc test/transport_parallel_stk.cpp)
        target_link_libraries(transport_parallel_stk ${amanzi_libs} ${UnitTest_LIBRARIES})
        add_amanzi_test(transport_parallel_stk transport_parallel_stk NPROCS 4 KIND uint)
    endif()    
endif()


