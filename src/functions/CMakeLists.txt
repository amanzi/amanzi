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
include_directories(${Tpetra_INCLUDE_DIRS})
include_directories(${HDF5_INCLUDE_DIRS})

# Need to move this up to the root -- lpritch
add_definitions("-DUSE_MPI")

#
# Library: functions
#
set(funcs_inc_files Function.hh
  FunctionComposition.hh
  FunctionMultiplicative.hh
  FunctionSmoothStep.hh
  FunctionStaticHead.hh
  FunctionAdditive.hh
  FunctionConstant.hh
  FunctionLinear.hh
  FunctionPolynomial.hh
  FunctionSquareDistance.hh
  FunctionTabular.hh
  FunctionBilinear.hh
  FunctionDistance.hh
  FunctionMonomial.hh
  FunctionSeparable.hh
  FunctionStandardMath.hh
  FunctionFactory.hh
  MultiFunction.hh
  )

set(funcs_src_files FunctionBilinear.cc
  FunctionMonomial.cc
  FunctionSmoothStep.cc
  FunctionStandardMath.cc
  FunctionDistance.cc
  FunctionLinear.cc
  FunctionPolynomial.cc
  FunctionSquareDistance.cc
  FunctionTabular.cc
  FunctionFactory.cc
  MultiFunction.cc
  )

set(funcs_link_libs error_handling
  ${HDF5_LIBRARIES}
  ${Teuchos_LIBRARIES}
  ${Tpetra_LIBRARIES}
  )

add_amanzi_library(functions 
                   SOURCE ${funcs_src_files}
                   HEADERS ${funcs_inc_files}
                   LINK_LIBS ${funcs_link_libs}
                   )

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

  add_amanzi_test(functions_all functions_all
                  KIND unit
                  SOURCE test/functions_all.cc
                  LINK_LIBS functions ${UnitTest_LIBRARIES})

  add_amanzi_test(functions_factory functions_factory
                  KIND unit
                  SOURCE test/functions_factory.cc
                  LINK_LIBS functions ${UnitTest_LIBRARIES})
 
  add_amanzi_test(functions_multi functions_multi
                  KIND unit
                  SOURCE test/functions_multi.cc
               	  LINK_LIBS functions ${UnitTest_LIBRARIES})

endif()
