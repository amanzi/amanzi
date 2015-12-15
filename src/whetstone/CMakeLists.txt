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
include_directories(${Epetra_INCLUDE_DIRS})

#
# Library: whetstone
#
add_amanzi_library(whetstone
                   SOURCE mfd3d.cc
                          mfd3d_diffusion.cc mfd3d_diffusion_surface.cc 
                          mfd3d_diffusion_other.cc 
                          mfd3d_elasticity.cc 
                          mfd3d_electromagnetics.cc 
                          nlfv.cc Tensor.cc DenseMatrix.cc
                   HEADERS DenseMatrix.hh
                           DenseVector.hh
                           lapack.hh
                           mfd3d.hh
                           mfd3d_diffusion.hh
                           mfd3d_elasticity.hh
                           mfd3d_electromagnetics.hh
                           nlfv.hh
                           Tensor.hh
                           WhetStoneDefs.hh
                           WhetStone_typedefs.hh
                   LINK_LIBS geometry mesh error_handling ${Epetra_LIBRARIES})


if ( BUILD_TESTS )
    # Add UnitTest includes
    include_directories(${UnitTest_INCLUDE_DIRS})
    include_directories(${MESH_SIMPLE_SOURCE_DIR})

    # Copy test directory files if an out of source build
    if (NOT (${WHETSTONE_SOURCE_DIR} EQUAL ${WHETSTONE_BINARY_DIR}) )
        execute_process(COMMAND ${CMAKE_COMMAND} -E 
          copy_directory ${WHETSTONE_SOURCE_DIR}/test ${WHETSTONE_BINARY_DIR}/test) 
    endif()

    set(amanzi_libs geometry mesh mesh_audit mesh_factory whetstone)

    # Test: matrices for diffusion
    add_amanzi_test(whetstone_darcy whetstone_darcy
                    KIND unit
                    SOURCE test/Main.cc test/whetstone_darcy.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    # Test: matrices for elasticity
    add_amanzi_test(whetstone_elasticity whetstone_elasticity
                    KIND unit
                    SOURCE test/Main.cc test/whetstone_elasticity.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    # Test: matrices for electromagnetics
    add_amanzi_test(whetstone_electromagnetics whetstone_electromagnetics
                    KIND unit
                    SOURCE test/Main.cc test/whetstone_electromagnetics.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    # Test: various operations with tensors 
    add_amanzi_test(whetstone_tensors whetstone_tensors
                    KIND unit
                    SOURCE test/Main.cc test/whetstone_tensors.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    # Test: surface matrices for diffusion
    add_amanzi_test(whetstone_surface whetstone_surface
                    KIND unit
                    SOURCE test/Main.cc test/whetstone_surface.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    # Test: various toolsets
    add_amanzi_test(whetstone_nlfv whetstone_nlfv
                    KIND unit
                    SOURCE test/Main.cc test/whetstone_nlfv.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})
endif()


