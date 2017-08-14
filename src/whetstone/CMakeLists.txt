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
include_directories(${Boost_INCLUDE_DIRS})

#
# Library: whetstone
#
add_amanzi_library(whetstone
                   SOURCE InnerProduct.cc 
                          DeRham_Node.cc DeRham_Edge.cc DeRham_Face.cc
                          MFD3D.cc
                          MFD3D_Diffusion.cc MFD3D_Diffusion_Surface.cc MFD3D_Diffusion_Other.cc 
                          MFD3D_Elasticity.cc 
                          MFD3D_Electromagnetics.cc MFD3D_Electromagnetics_Surface.cc
                          MFD3D_Generalized_Diffusion.cc
                          MFD3D_BernardiRaugel.cc
                          DG_Modal.cc nlfv.cc 
                          Tensor.cc DenseMatrix.cc
                          Polynomial.cc
                          MeshMaps.cc
                          MeshMaps_FEM.cc MeshMaps_VEM.cc
                   HEADERS DenseMatrix.hh DenseVector.hh
                           DG_Modal.hh
                           lapack.hh
                           InnerProductL2.hh InnerProductH1.hh BilinearForm.hh
                           DeRham_Node.hh DeRham_Edge.hh DeRham_Face.hh
                           MFD3D.hh
                           MFD3D_Diffusion.hh
                           MFD3D_Elasticity.hh
                           MFD3D_Electromagnetics.hh
                           MFD3D_Generalized_Diffusion.hh
                           MFD3D_BernardiRaugel.hh
                           MFD3DFactory.hh
                           nlfv.hh
                           Polynomial.hh Tensor.hh
                           MeshMaps.hh
                           MeshMaps_FEM.hh MeshMaps_VEM.hh
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
    add_amanzi_test(whetstone_diffusion whetstone_diffusion
                    KIND unit
                    SOURCE test/Main.cc test/whetstone_diffusion.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    # Test: matrices for elasticity
    add_amanzi_test(whetstone_elasticity whetstone_elasticity
                    KIND unit
                    SOURCE test/Main.cc test/whetstone_elasticity.cc
                           test/whetstone_elasticity_other.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    # Test: matrices for electromagnetics
    add_amanzi_test(whetstone_electromagnetics whetstone_electromagnetics
                    KIND unit
                    SOURCE test/Main.cc test/whetstone_electromagnetics.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    # Test: various operations with tensors 
    add_amanzi_test(whetstone_tensors_geometry whetstone_tensors_geometry
                    KIND unit
                    SOURCE test/Main.cc test/whetstone_tensors.cc
                                        test/whetstone_geometry.cc
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

    # Test: matrices for DG methods
    add_amanzi_test(whetstone_dg whetstone_dg
                    KIND unit
                    SOURCE test/Main.cc test/whetstone_dg.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    # Test: matrices for Derham complex
    add_amanzi_test(whetstone_derham whetstone_derham
                    KIND unit
                    SOURCE test/Main.cc test/whetstone_derham.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

endif()


