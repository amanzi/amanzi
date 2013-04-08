# -*- mode: cmake -*-

#
#  Amanzi
#    Mesh functions
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(LibraryManager)

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
include_directories(${DATA_STRUCTURES_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})

#
# Library: mesh_functions
#
file(GLOB src_files "*.cc")
file(GLOB inc_files "*.hh")
add_amanzi_library(mesh_functions
                   SOURCE ${src_files} HEADERS ${inc_files}
                   LINK_LIBS error_handling functions mesh ${Teuchos_LIBRARIES})

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

# these are currently broken
#  add_amanzi_test(boundary_function_test boundary_function_test
#                  KIND unit
#                  SOURCE test/boundary-function-test.cc
#                  LINK_LIBS mesh_functions mesh_factory ${UnitTest_LIBRARIES})

  add_amanzi_test(cv_function cv_function_test 
                  KIND unit
                  SOURCE test/composite_vector_function_test.cc
                  LINK_LIBS  mesh_functions mesh_factory data_structures 
                             ${UnitTest_LIBRARIES} ${Boost_LIBRARIES} 
                             ${MSTK_LIBRARIES})



endif()
