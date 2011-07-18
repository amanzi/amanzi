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

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${NOX_INCLUDE_DIRS})

# Need to move this to roo level -- lpritch
add_definitions("-DUSE_MPI")

#
# Library: time_integration
#
file(GLOB time_integration_source_files "*.cpp")
add_library(time_integration ${time_integration_source_files})
add_install_library(time_integration)
link_list_add(time_integration)
set(time_integration_tpl_libs ${DBC_LIBRARIES} ${Teuchos_LIBRARIES} ${Epetra_LIBRARIES} ${NOX_LIBRARIES})
list(REMOVE_DUPLICATES time_integration_tpl_libs)
target_link_libraries(time_integration error_handling flow ${time_integration_tpl_libs})

#
# Install Targets
#
file(GLOB time_int_inc_files "*.hpp")
add_install_include_file(${time_int_inc_files})


if (BUILD_TESTS)

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

    add_executable(test_time_integration test/Main.cc test/test_solhist.cpp test/test_1Dfem.cpp) 
    target_link_libraries(test_time_integration time_integration solvers error_handling ${UnitTest_LIBRARIES})
    add_amanzi_test(time_integration test_time_integration KIND unit)

endif()
