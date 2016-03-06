# -*- mode: cmake -*-

#
#  Amanzi
#   Operator library
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)
include(LibraryManager)

#
# Define a project name
# After this command the following varaibles are defined
#   OPERATORS_SOURCE_DIR
#   OPERATORS_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(OPERATORS)

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${MSTK_INCLUDE_DIRS})
# Amanzi include directories
include_directories(${ATK_SOURCE_DIR})
include_directories(${DBC_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${MESH_LOGICAL_SOURCE_DIR})
include_directories(${MSTK_SOURCE_DIR})
include_directories(${MESH_FACTORY_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${WHETSTONE_SOURCE_DIR})
include_directories(${DATA_STRUCTURES_SOURCE_DIR})
include_directories(${SOLVERS_SOURCE_DIR})
include_directories(${SOLVERS_SOURCE_DIR})
include_directories(${OPERATORS_SOURCE_DIR})
include_directories(${OUTPUT_SOURCE_DIR})

# Need this define. Errors from MSTK include files 
# about MPI_COMM_WORLD. --lpritch
add_definitions("-DMSTK_HAVE_MPI")

#
# Library: operator audit
#
add_amanzi_library(operators_audit SOURCE OperatorAudit.cc LINK_LIBS ${flow_tpl_libs})

#
# Library: operators
#
set(operators_inc_files
    BCs.hh
    FieldMaps.hh
    OperatorDefs.hh
    OperatorUtils.hh
    Op.hh
    Op_Cell_Cell.hh
    Op_Cell_Face.hh
    Op_Cell_FaceCell.hh
    Op_Cell_Node.hh
    Op_Cell_Edge.hh
    Op_Face_Cell.hh
    Op_Node_Node.hh
    Operator.hh
    Operator_Cell.hh
    Operator_FaceCell.hh
    Operator_FaceCellScc.hh
    Operator_FaceCellSff.hh
    Operator_ConsistentFace.hh
    Operator_Node.hh
    Operator_Edge.hh
    BoundaryFlux.hh
    OperatorCurlCurl.hh
    OperatorDiffusion.hh
    OperatorDiffusionWithGravity.hh
    OperatorDiffusionFactory.hh
    OperatorDiffusionMFD.hh
    OperatorDiffusionMFDwithGravity.hh
    OperatorDiffusionFV.hh
    OperatorDiffusionFVwithGravity.hh
    OperatorDiffusionNLFV.hh
    OperatorDiffusionNLFVwithGravity.hh
    OperatorAccumulation.hh
    OperatorAdvection.hh
    TreeOperator.hh
    Reconstruction.hh
    ReconstructionCell.hh
    Upwind.hh
    UpwindArithmeticAverage.hh
    UpwindDivK.hh
    UpwindFactory.hh
    UpwindFlux.hh
    UpwindGravity.hh
    UpwindSecondOrder.hh
)

add_amanzi_library(operators 
                   SOURCE OperatorUtils.cc
                          FieldMaps.cc
                          Operator.cc
                          Operator_Cell.cc Operator_FaceCell.cc
                          Operator_FaceCellScc.cc
                          Operator_FaceCellSff.cc
                          Operator_ConsistentFace.cc
                          Operator_Node.cc
                          Operator_Edge.cc
                          OperatorCurlCurl.cc
                          OperatorDiffusionMFD.cc
                          OperatorDiffusionMFDwithGravity.cc
                          OperatorDiffusionFV.cc
                          OperatorDiffusionFVwithGravity.cc
                          OperatorDiffusionNLFV.cc
                          OperatorDiffusionNLFVwithGravity.cc
                          OperatorDiffusionFactory.cc                   
                          OperatorAccumulation.cc
                          OperatorAdvection.cc  
                          TreeOperator.cc
                          ReconstructionCell.cc
                          ReconstructionCell_Limiters.cc
                   HEADERS ${operators_inc_files}
                   LINK_LIBS geometry mesh error_handling whetstone solvers data_structures ${Epetra_LIBRARIES})


if ( BUILD_TESTS )
    # Add UnitTest includes
    include_directories(${UnitTest_INCLUDE_DIRS})
    include_directories(${MESH_SIMPLE_SOURCE_DIR})

    # Copy test directory files if an out of source build
    if (NOT (${OPERATORS_SOURCE_DIR} EQUAL ${OPERATORS_BINARY_DIR}) )
        file(GLOB DataFiles "${OPERATORS_SOURCE_DIR}/test/*.xml"
                            "${OPERATORS_SOURCE_DIR}/test/*.exo")
        file(COPY ${DataFiles} DESTINATION ${OPERATORS_BINARY_DIR}/test/)
    endif()

    set(amanzi_libs geometry atk mesh mesh_factory mesh_logical data_structures operators whetstone output atk)

    # Test: operators for various PDEs
    if ( ENABLE_INDIVIDUAL_TESTS_OPERATORS OR ENABLE_INDIVIDUAL_TESTS )
        add_amanzi_test(operators_laplace_beltrami operators_laplace_beltrami
                        KIND unit
                        SOURCE test/Main.cc test/operator_laplace_beltrami.cc 
                        LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})
        add_amanzi_test(operators_parabolic_surface operators_parabolic_surface
                        KIND unit
                        SOURCE test/Main.cc test/operator_parabolic_surface.cc
                        LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})
        add_amanzi_test(operators_marshak operators_marshak
                        KIND unit
                        SOURCE test/Main.cc test/operator_marshak.cc
                        LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})
        add_amanzi_test(operators_marshak_logical operators_marshak_logical
                        KIND unit
                        SOURCE test/Main.cc test/operator_marshak_logical.cc
                        LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})
    else()
        add_amanzi_test(operators_bundle_pdes operators_bundle_pdes
                        KIND unit
                        SOURCE test/Main.cc test/operator_laplace_beltrami.cc 
                                            test/operator_parabolic_surface.cc
                                            test/operator_marshak.cc
                                            test/operator_marshak_logical.cc
                        LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})
    endif()

    add_amanzi_test(operators_nonlinear_surface operators_nonlinear_surface
                    KIND unit
                    SOURCE test/Main.cc test/operator_nonlinear_surface.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    add_amanzi_test(operators_diffusion operators_diffusion
                    KIND unit
                    SOURCE test/Main.cc test/operator_diffusion.cc
                                        test/operator_diffusion_nodal.cc
                                        test/operator_diffusion_divk.cc
                                        test/operator_diffusion_nonsymmetric.cc
                                        test/operator_diffusion_nlfv.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    add_amanzi_test(operators_diffusion_3pe operators_diffusion NPROCS 3 KIND unit)

    add_amanzi_test(operators_diffusion_consistent_face operators_diffusion_consistent_face
                    KIND unit
                    SOURCE test/Main.cc test/operator_diffusion_consistent_face.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    add_amanzi_test(operators_diffusion_consistent_face_3pe operators_diffusion_consistent_face NPROCS 3 KIND unit)

    add_amanzi_test(operators_commute operators_commute
                    KIND unit
                    SOURCE test/Main.cc test/operator_commute.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    # # Test: tpfa
    # add_amanzi_test(operators_diffusion_tpfa operators_diffusion_tpfa
    #                 KIND unit
    #                 SOURCE test/Main.cc test/operator_diffusion_tpfa.cc
    #                 LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})
    
    # Test: matrices for surface
    add_amanzi_test(operators_advdiff_surface operators_advdiff_surface
                    KIND unit
                    SOURCE test/Main.cc test/operator_advdiff_surface.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    # Test: operators on 
    add_amanzi_test(operators_mhd operators_mhd
                    KIND unit
                    SOURCE test/Main.cc test/operator_mhd.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    # Test: utils
    add_amanzi_test(operators_upwind operators_upwind
                    KIND unit
                    SOURCE test/Main.cc test/operator_upwind.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    add_amanzi_test(operators_upwind_2pe operators_upwind NPROCS 2 KIND unit)

    # stability analysis
    add_amanzi_test(operators_stability operators_stability
                    KIND unit
                    SOURCE test/Main.cc test/operator_stability.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    add_amanzi_test(operators_stability_2pe operators_stability NPROCS 2 KIND unit)

    # correctness of tree operator
    add_amanzi_test(operators_tree operators_tree
                    KIND unit
                    SOURCE test/Main.cc test/operator_tree.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    add_amanzi_test(operators_tree_2pe operators_tree NPROCS 2 KIND unit)

    # Test: diffusion with gravity
    add_amanzi_test(operators_gravity operators_gravity
                    KIND unit
                    SOURCE test/Main.cc test/operator_gravity.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    # reconstruction tests
    add_amanzi_test(operators_reconstruction operators_reconstruction
                    KIND unit
                    SOURCE test/Main.cc test/operator_reconstruction.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    add_amanzi_test(operators_boundary_flux_solver operators_boundary_flux_solver
                    KIND unit
                    SOURCE test/Main.cc test/operator_boundary_flux_solver.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    add_amanzi_test(operators_reconstruction_2pe operators_reconstruction NPROCS 2 KIND unit)

endif()


