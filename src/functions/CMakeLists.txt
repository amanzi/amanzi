# -*- mode: cmake -*-

#
#  Amanzi
#    Functions
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)
include(LibraryManager)
#include(RegistryManager)

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
include_directories(${ATK_SOURCE_DIR})
include_directories(${DBC_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${HDF5_INCLUDE_DIRS})

# Need to move this up to the root -- lpritch
add_definitions("-DUSE_MPI")

#
# Library: functions
#
file(GLOB src_files "*.cc")
file(GLOB inc_files "*.hh")
add_amanzi_library(functions 
                   SOURCE ${src_files}
                   HEADERS ${inc_files}
                   LINK_LIBS atk error_handling ${HDF5_LIBRARIES} ${Epetra_LIBRARIES} ${Teuchos_LIBRARIES})

# Define all tests at the of the file. BUILD_TESTS gates test building.
if (BUILD_TESTS)

  # Add UnitTest includes
  include_directories(${UnitTest_INCLUDE_DIRS})

  # Add include directory. Need to remove this required path. -- lpritch
  include_directories(${FUNCS_SOURCE_DIR})

  # Copy test subdirectory for out of source builds
  if (NOT (FUNCS_SOURCE_DIR STREQUAL FUNCS_BINARY_DIR))
      file(GLOB DataFiles "${FUNCS_SOURCE_DIR}/test/*.txt"
                          "${FUNCS_SOURCE_DIR}/test/*.h5")
      file(COPY ${DataFiles} DESTINATION ${FUNCS_BINARY_DIR}/test/)
  endif()

  add_amanzi_test(function_all function_all
                  KIND unit
                  SOURCE test/AllFunctionTest.cc
                  LINK_LIBS functions ${UnitTest_LIBRARIES})

  add_amanzi_test(function_factory function_factory
                  KIND unit
                  SOURCE test/FunctionFactoryTest.cc
                  LINK_LIBS functions ${UnitTest_LIBRARIES})
 
  add_amanzi_test(function_concrete_color function_concrete_color
                  KIND unit
                  SOURCE test/AllColorFunctionTest.cc
                  LINK_LIBS functions ${UnitTest_LIBRARIES})

  add_amanzi_test(function_color_factory function_color_factory
                  KIND unit
                  SOURCE test/ColorFunctionFactoryTest.cc
                  LINK_LIBS functions ${UnitTest_LIBRARIES})

  add_amanzi_test(function_composite function_composite
                  KIND unit
                  SOURCE test/CompositeFunctionTest.cc
		  LINK_LIBS functions ${UnitTest_LIBRARIES})

endif()
