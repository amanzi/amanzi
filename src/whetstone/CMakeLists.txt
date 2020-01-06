# -*- mode: cmake -*-

#
#  Amanzi
#   Discretization library
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)
include(LibraryManager)

#
# Define a project name
# After this command the following varaibles are defined
#   WHETSTONE_SOURCE_DIR
#   WHETSTONE_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(WHETSTONE)

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${MESH_FACTORY_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${WHETSTONE_SOURCE_DIR})
include_directories(${DBG_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Tpetra_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})

set(whetstone_inc_files
  lapack.hh
  WhetStoneDefs.hh
  WhetStoneFunction.hh
  WhetStoneMeshUtils.hh

  DenseVector.hh
  DenseMatrix.hh
  Tensor.hh

  PolynomialBase.hh
  PolynomialIterator.hh
  Monomial.hh
  Polynomial.hh
  SpaceTimePolynomial.hh
  VectorObjects.hh
  MatrixObjects.hh
)

set(whetstone_src_files
  DenseVector.cc
  DenseMatrix.cc
  Tensor.cc
  Polynomial.cc
  SpaceTimePolynomial.cc
  VectorObjects.cc
)

#
# Library: whetstone
#
add_amanzi_library(whetstone
  SOURCE ${whetstone_src_files}
  HEADERS ${whetstone_inc_files}
  LINK_LIBS ${whetstone_link_libs}
  )


if (BUILD_TESTS)
    # Add UnitTest includes
    include_directories(${UnitTest_INCLUDE_DIRS})

    # Copy test directory files if an out of source build
    if (NOT (${WHETSTONE_SOURCE_DIR} EQUAL ${WHETSTONE_BINARY_DIR}))
        file(GLOB whetstone_data_files "${WHETSTONE_SOURCE_DIR}/test/*.exo")
        file(COPY ${whetstone_data_files} DESTINATION ${WHETSTONE_BINARY_DIR}/test/)
    endif()

    set(whetstone_test_link_libs geometry whetstone)

    # Test: miscalleneous library tools
    add_amanzi_test(whetstone_misc whetstone_misc
                    KIND unit
                    SOURCE test/Main.cc test/whetstone_tensors.cc
                    LINK_LIBS ${whetstone_test_link_libs} ${UnitTest_LIBRARIES})

endif()


