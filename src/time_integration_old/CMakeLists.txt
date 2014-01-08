# -*- mode: cmake -*-

#
#  Amanzi
#   Flow process kernel
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)

#
# Define a project name
# After this command the following varaibles are defined
#   TIME_INTEGRATION_SOURCE_DIR
#   TIME_INTEGRATION_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(TIME_INTEGRATION)

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})
include_directories(${SOLVERS_SOURCE_DIR})
include_directories(${FLOW_SOURCE_DIR})
include_directories(${WHETSTONE_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${STATE_SOURCE_DIR})
include_directories(${FUNCS_SOURCE_DIR})
include_directories(${MFUNCS_SOURCE_DIR})
include_directories(${OUTPUT_SOURCE_DIR})
include_directories(${PLATFORM_INTERFACE_SOURCE_DIR})
include_directories(${DBG_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${NOX_INCLUDE_DIRS})
include_directories(${HDF5_C_INCLUDE_DIR})
include_directories(${ASCEMIO_INCLUDE_DIR})

#
# Library: time_integration
#
file(GLOB time_integration_source_files "*.cc")
file(GLOB time_int_inc_files "*.hh")
set(time_integration_tpl_libs ${DBC_LIBRARIES} ${Teuchos_LIBRARIES} ${Epetra_LIBRARIES} ${NOX_LIBRARIES})
list(REVERSE time_integration_tpl_libs)
list(REMOVE_DUPLICATES time_integration_tpl_libs)
list(REVERSE time_integration_tpl_libs)
add_amanzi_library(time_integration
                   SOURCE ${time_integration_source_files} HEADERS ${time_int_inc_files}
                   LINK_LIBS error_handling solvers ${time_integration_tpl_libs} atk)
if (BUILD_TESTS)
#if (0)

    # Add UnitTest include directoy
    include_directories(${UnitTest_INCLUDE_DIRS})

    # Copy test directory files if an out of source build
    if (NOT (${TIME_INTEGRATION_SOURCE_DIR} EQUAL ${TIME_INTEGRATION_BINARY_DIR}))
        execute_process(COMMAND ${CMAKE_COMMAND} -E 
          copy_directory ${TIME_INTEGRATION_SOURCE_DIR}/test ${TIME_INTEGRATION_BINARY_DIR}/test) 
    endif()

    # Add the time_integration directory to the include paths
    include_directories(${TIME_INTEGRATION_SOURCE_DIR})
    include_directories(${SOLVERS_SOURCE_DIR})

    add_amanzi_test(time_integration test_time_integration
                    KIND unit
                    SOURCE 
   		          test/Main.cc
   			  test/test_solhist.cc
   			  test/test_Explicit_RK.cc 
                    LINK_LIBS
   		          time_integration
   			  solvers
   			  error_handling
   			  ${UnitTest_LIBRARIES}
   			  ${time_integration_tpl_libs})

endif()
