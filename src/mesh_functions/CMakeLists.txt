# -*- mode: cmake -*-

#
#  Amanzi
#    Mesh functions
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
#include(AddParallelTest)

#
# Define a project name
# After this command the following variables are defined
#   XXX_SOURCE_DIR
#   XXX_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(MFUNCS)

# Define include directories to build any binary or library
# in this directory. 

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${FUNCS_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})

# Need to move this up to the root -- lpritch
add_definitions("-DUSE_MPI")

#
# Library: mesh_functions
#
file(GLOB src_files "*.cc")
add_library(mesh_functions ${src_files})
add_install_library(mesh_functions)
target_link_libraries(mesh_functions error_handling functions mesh ${Teuchos_LIBRARIES})

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
  include_directories(${MFUNCS_SOURCE_DIR})
  include_directories(${FUNCS_SOURCE_DIR})
  include_directories(${MESH_FACTORY_SOURCE_DIR})

  # Copy test subdirectory for out of source builds
  if (NOT (MFUNCS_SOURCE_DIR STREQUAL MFUNCS_BINARY_DIR))
      execute_process(COMMAND ${CMAKE_COMMAND} -E 
        copy_directory ${MFUNCS_SOURCE_DIR}/test ${MFUNCS_BINARY_DIR}/test) 
  endif()

  add_executable(boundary-function-test test/boundary-function-test.cc)
  target_link_libraries(boundary-function-test mesh_functions mesh_factory ${UnitTest_LIBRARIES})
  add_amanzi_test(boundary-function boundary-function-test KIND unit)

endif()
