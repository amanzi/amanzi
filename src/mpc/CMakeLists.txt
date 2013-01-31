# -*- mode: cmake -*-

#
#  Amanzi
#    Mulitprocess Coordinator (MPC)
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)
include(LibraryManager)

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
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${OUTPUT_SOURCE_DIR})
include_directories(${STATE_SOURCE_DIR})
include_directories(${CHEM_SOURCE_DIR})
include_directories(${CHEMPK_SOURCE_DIR})
include_directories(${TRANSPORT_SOURCE_DIR})
include_directories(${WHETSTONE_SOURCE_DIR})
include_directories(${FLOW_SOURCE_DIR})
include_directories(${TIME_INTEGRATION_SOURCE_DIR})
include_directories(${SOLVERS_SOURCE_DIR})
include_directories(${PLATFORM_INTERFACE_SOURCE_DIR})
include_directories(${FUNCS_SOURCE_DIR})
include_directories(${MFUNCS_SOURCE_DIR})
include_directories(${DBG_SOURCE_DIR})


# Chemistry includes directory
get_property(CHEM_INCLUDES_DIR GLOBAL PROPERTY CHEM_INCLUDES_DIR)
include_directories(${CHEM_INCLUDES_DIR})

# External (TPL) include directories
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${HDF5_C_INCLUDE_DIR})
include_directories(${ASCEMIO_INCLUDE_DIR})

#
# Library: mpc
#
set(boost_libs_used ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${Boost_REGEX_LIBRARY})
set(mpc_tpl_libs ${boost_libs_used} ${Teuchos_LIBRARIES} ${Epetra_LIBRARIES} ) 
list(REVERSE mpc_tpl_libs)
list(REMOVE_DUPLICATES mpc_tpl_libs) # Concerned about command line length! --lpritch
list(REVERSE mpc_tpl_libs)
file(GLOB mpc_inc_files "*.h")
add_amanzi_library(mpc 
                   SOURCE MPC.cc Unstructured_observations.cc HEADERS ${mpc_inc_files}
                   LINK_LIBS time_integration solvers output state mesh dbg ${mpc_tpl_libs})

if (BUILD_TESTS)
#if (0)

  include_directories(${MESH_FACTORY_SOURCE_DIR})

    
    # Add UnitTest includes
    include_directories(${UnitTest_INCLUDE_DIRS})

    # Add MPC include directory. Need to remove this required path. -- lpritch
    include_directories(${MPC_SOURCE_DIR})

    # Copy test subdirectory for out of source builds
    if (NOT (MPC_SOURCE_DIR STREQUAL MPC_BINARY_DIR))
        execute_process(COMMAND ${CMAKE_COMMAND} -E 
          copy_directory ${MPC_SOURCE_DIR}/test ${MPC_BINARY_DIR}/test) 
    endif()

    set(amanzi_libs geometry mesh mesh_audit mesh_factory state mpc flow transport chemistry whetstone)
    
    if ( ENABLE_MSTK_Mesh )
        # Test: mpc-mstk
	add_amanzi_test(mpc-mstk test_mpc_mstk
	                KIND      int
			SOURCE    test/Main.cc test/test_driver.cpp
			LINK_LIBS mpc
			          mesh
				  mesh_factory
				  chemistry_pk
				  flow
				  transport
				  whetstone
				  ${UnitTest_LIBRARIES})

    
    endif()

#    add_executable(test_vis test/Main.cc test/test_vis.cpp)
#    target_link_libraries(test_vis state mpc stk_mesh ${UnitTest_LIBRARIES})
#    add_amanzi_test(vis test_vis KIND int)

endif()
