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
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${ASCEMIO_INCLUDE_DIR})


# Need to move this to the root level -- lpritch
add_definitions("-DUSE_MPI" "-DMSTK_HAVE_MPI")

#
# Library: transport
#
file(GLOB transport_inc_files "*.hpp")
set(transport_src_files Transport_PK.cpp Transport_State.cpp Transport_TI.cpp 
                        Transport_Limiters.cpp Transport_VandV.cpp Transport_IO.cpp
                        Reconstruction.cpp)
add_amanzi_library(transport 
                   SOURCE ${transport_src_files} HEADERS ${transport_inc_files}
                   LINK_LIBS geometry mesh mesh_data error_handling state time_integration ${Epetra_LIBRARIES})

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

    set(amanzi_libs geometry mesh mesh_audit simple_mesh state transport whetstone mesh_functions)

    # Test: miscaleneous transport routines
    add_amanzi_test(transport_mics transport_mics
                    KIND unit
                    SOURCE test/Main.cc test/transport_mics.cpp
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    add_amanzi_test(transport_advance_simple transport_advance_simple
                    KIND unit
                    SOURCE test/Main.cc test/transport_advance_simple.cpp
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    add_amanzi_test(transport_convergence transport_convergence
                    KIND unit
                    SOURCE test/Main.cc test/transport_convergence.cpp
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    # add_executable(transport_dispersion test/Main.cc test/transport_dispersion.cpp)
    # target_link_libraries(transport_dispersion ${amanzi_libs} ${UnitTest_LIBRARIES})
    # add_amanzi_test(transport_dispersion transport_dispersion KIND unit)

    if ( ENABLE_MSTK_Mesh )
        include_directories(${MSTK_INCLUDE_DIRS})
        include_directories(${MSTK_SOURCE_DIR})

        set(amanzi_libs geometry mesh mesh_audit mstk_mesh state transport whetstone output mesh_functions)

        add_amanzi_test(transport_advance_mstk transport_advance_mstk
	                KIND unit
                        SOURCE test/Main.cc test/transport_advance_mstk.cpp
                        LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

        add_amanzi_test(transport_parallel_mstk transport_parallel_mstk
		       KIND int
	               NPROCS 4 
                       SOURCE test/Main.cc test/transport_parallel_mstk.cpp
                       LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

        add_amanzi_test(transport_2D transport_2D
	                KIND unit
                        SOURCE test/Main.cc test/transport_2D.cpp
                        LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

        add_amanzi_test(transport_2D_long transport_2D_long 
	                KIND unit
                        SOURCE test/Main.cc test/transport_2D_long.cpp
                        LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

        add_amanzi_test(transport_3D_long transport_3D_long
	                KIND unit
                        SOURCE test/Main.cc test/transport_3D_long.cpp
                        LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

        add_amanzi_test(transport_subcycling transport_subcycling
	                KIND unit
                        SOURCE test/Main.cc test/transport_subcycling.cpp
                        LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})
    endif()    

    if ( ENABLE_STK_Mesh )
        include_directories(${STK_INCLUDE_DIRS})
        include_directories(${STK_SOURCE_DIR})

        set(amanzi_libs geometry mesh mesh_audit stk_mesh exodus_reader state transport whetstone mesh_functions)

        add_amanzi_test(transport_advance_stk transport_advance_stk
	                KIND unit
                        SOURCE test/Main.cc test/transport_advance_stk.cpp
                        LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})
        
        add_amanzi_test(transport_parallel_stk transport_parallel_stk
			KIND unit
	                NPROCS 4
                        SOURCE test/Main.cc test/transport_parallel_stk.cpp
                        LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})
          
    endif()    
endif()


