# -*- mode: cmake -*-

#
#  Amanzi
#    Functions
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)

#
# Define a project name
# After this command the following variables are defined
#   XXX_SOURCE_DIR
#   XXX_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(FUNCS)

# Define include directories to build any binary or library
# in this directory. 

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})

# Need to move this up to the root -- lpritch
add_definitions("-DUSE_MPI")

#
# Library: functions
#
file(GLOB src_files "*.cc")
add_library(functions ${src_files})
add_install_library(functions)
target_link_libraries(functions error_handling ${Teuchos_LIBRARIES} ${Epetra_LIBRARIES})

#
# Install Targets
#
file(GLOB inc_files "*.hh")
add_install_include_file(${inc_files})

# Define all tests at the of the file. BUILD_TESTS gates test building.
if (BUILD_TESTS)
    
  # Add UnitTest includes
  include_directories(${UnitTest_INCLUDE_DIRS})

  # Add include directory. Need to remove this required path. -- lpritch
  include_directories(${FUNCS_SOURCE_DIR})

  # Copy test subdirectory for out of source builds
  if (NOT (FUNCS_SOURCE_DIR STREQUAL FUNCS_BINARY_DIR))
      execute_process(COMMAND ${CMAKE_COMMAND} -E 
        copy_directory ${FUNCS_SOURCE_DIR}/test ${FUNCS_BINARY_DIR}/test) 
  endif()

  add_amanzi_test(concrete-functions all-function-test
                  KIND unit
  SOURCE test/all-function-test.cc
  LINK_LIBS functions ${UnitTest_LIBRARIES})

  add_amanzi_test(function-factory function-factory-test
                  KIND unit
                  SOURCE test/function-factory-test.cc
                  LINK_LIBS functions ${UnitTest_LIBRARIES})
 
  add_amanzi_test(concrete-color-function all-color-function-test
                  KIND unit
                  SOURCE test/all-color-function-test.cc
                  LINK_LIBS functions ${UnitTest_LIBRARIES})

  add_amanzi_test(color-function-factory color-function-factory-test
                  KIND unit
                  SOURCE test/color-function-factory-test.cc
                  LINK_LIBS 
                         functions
                         ${UnitTest_LIBRARIES}
                         ${Teuchos_LIBRARIES}
                         ${Epetra_LIBRARIES})

endif()
