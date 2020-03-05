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
include_directories(${Tpetra_INCLUDE_DIRS})
include_directories(${Trilinos_INCLUDE_DIRS})

# Amanzi include directories
include_directories(${ATK_SOURCE_DIR})
include_directories(${DATA_STRUCTURES_SOURCE_DIR})
include_directories(${DBC_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
#include_directories(${MESH_LOGICAL_SOURCE_DIR})
include_directories(${MESH_FACTORY_SOURCE_DIR})
include_directories(${MSTK_SOURCE_DIR})
include_directories(${OPERATORS_SOURCE_DIR})
include_directories(${OUTPUT_SOURCE_DIR})
include_directories(${SOLVERS_SOURCE_DIR})
include_directories(${TIME_INTEGRATION_SOURCE_DIR})
include_directories(${WHETSTONE_SOURCE_DIR})


# Need this define. Errors from MSTK include files 
# about MPI_COMM_WORLD. --lpritch
add_definitions("-DMSTK_HAVE_MPI")

#
# Library: operators
#
set(operators_inc_files
  OperatorDefs.hh
  Schema.hh
  Op.hh
  Op_Cell_Cell.hh
  Op_Face_Cell.hh
  Op_Face_CellBndFace.hh
  Op_SurfaceCell_SurfaceCell.hh  
  Op_SurfaceFace_SurfaceCell.hh
  OperatorUtils.hh
  Operator.hh
  Operator_Cell.hh
  Operator_CellBndFace.hh
  Operator_Factory.hh
  TreeOperator.hh
  BCs.hh
  PDE_HelperBCsList.hh
  PDE_HelperDiscretization.hh
  PDE_Diffusion.hh
  PDE_DiffusionFV.hh
  PDE_DiffusionFVwithGravity.hh
  )

set(operators_src_files
  Schema.cc
  Op.cc
  OperatorUtils.cc
  Operator.cc
  Operator_Cell.cc
  Operator_CellBndFace.cc
  Operator_Factory.cc
  TreeOperator.cc
  PDE_HelperDiscretization.cc
  PDE_DiffusionFV.cc
  PDE_DiffusionFVwithGravity.cc
  )

set(operators_link_libs
  error_handling
  whetstone
  mesh
  data_structures
  solvers
  ${Teuchos_LIBRARIES}
  ${Tpetra_LIBRARIES}
  )

# geometry  mesh_factory  solvers 

#   BCs.hh
#     LimiterCell.hh
#     Mini_Diffusion1D.hh
#     Mini_Operator1D.hh
#     OperatorDefs.hh
#     OperatorUtils.hh
#     Op.hh
#     Op_Cell_Cell.hh
#     Op_Cell_Face.hh
#     Op_Cell_FaceCell.hh
#     Op_Cell_Node.hh
#     Op_Cell_Edge.hh
#     Op_Cell_Schema.hh
#     Op_Diagonal.hh
#     Op_Face_Cell.hh
#     Op_Face_CellBndFace.hh
#     Op_Face_Schema.hh
#     Op_Node_Node.hh
#     Operator.hh
#     Operator_Cell.hh
#     Operator_CellBndFace.hh
#     Operator_ConsistentFace.hh
#     Operator_Diagonal.hh
#     Operator_Edge.hh
#     Operator_FaceCell.hh
#     Operator_FaceCellScc.hh
#     Operator_FaceCellSff.hh
#     Operator_Factory.hh
#     Operator_Node.hh
#     BoundaryFlux.hh
#     PDE_Abstract.hh
#     PDE_Accumulation.hh
#     PDE_Advection.hh
#     PDE_AdvectionUpwind.hh
#     PDE_AdvectionRiemann.hh
#     PDE_CouplingFlux.hh
#     PDE_Diffusion.hh
#     PDE_DiffusionWithGravity.hh
#     PDE_DiffusionFactory.hh
#     PDE_DiffusionMFD.hh
#     PDE_DiffusionMFDwithGravity.hh
#     PDE_DiffusionFV.hh
#     PDE_DiffusionFVwithGravity.hh
#     PDE_DiffusionNLFV.hh
#     PDE_DiffusionNLFVwithBndFaces.hh
#     PDE_DiffusionNLFVwithGravity.hh
#     PDE_DiffusionDG.hh
#     PDE_DiffusionFracturedMatrix.hh
#     PDE_Elasticity.hh
#     PDE_Electromagnetics.hh
#     PDE_MagneticDiffusion.hh
#     PDE_MagneticDiffusion_TM.hh
#     PDE_HelperBCsList.hh
#     PDE_HelperDiscretization.hh
#     PDE_Reaction.hh
#     RemapDG.hh
#     RemapUtils.hh
#     Schema.hh
#     TreeOperator.hh
#     Reconstruction.hh
#     ReconstructionCell.hh
#     Upwind.hh
#     UpwindArithmeticAverage.hh
#     UpwindDivK.hh
#     UpwindFactory.hh
#     UpwindFlux.hh
#     UpwindGravity.hh
#     UpwindSecondOrder.hh
# )



# LimiterCell.cc
#                           Mini_Diffusion1D.cc
#                           Mini_Operator1D.cc
#                           OperatorUtils.cc
#                           Operator.cc
#                           Operator_ConsistentFace.cc
#                           Operator_Edge.cc
#                           Operator_Cell.cc Operator_CellBndFace.cc
#                           Operator_Diagonal.cc
#                           Operator_FaceCell.cc
#                           Operator_FaceCellScc.cc Operator_FaceCellSff.cc
#                           Operator_Factory.cc
#                           Operator_Node.cc
#                           Operator_Schema.cc
#                           PDE_Accumulation.cc
#                           PDE_Abstract.cc
#                           PDE_AdvectionUpwind.cc PDE_AdvectionRiemann.cc
#                           PDE_CouplingFlux.cc
#                           PDE_DiffusionMFD.cc PDE_DiffusionMFDwithGravity.cc
#                           PDE_DiffusionFV.cc PDE_DiffusionFVwithGravity.cc
#                           PDE_DiffusionNLFV.cc PDE_DiffusionNLFVwithGravity.cc
#                           PDE_DiffusionNLFVwithBndFaces.cc PDE_DiffusionNLFVwithBndFacesGravity.cc
#                           PDE_DiffusionDG.cc
#                           PDE_DiffusionFracturedMatrix.cc
#                           PDE_DiffusionFactory.cc                   
#                           PDE_HelperDiscretization.cc
#                           PDE_Elasticity.cc
#                           PDE_Electromagnetics.cc
#                           PDE_MagneticDiffusion.cc PDE_MagneticDiffusion_TM.cc
#                           PDE_Reaction.cc
#                           Schema.cc
#                           TreeOperator.cc
#                           ReconstructionCell.cc
#                           RemapDG.cc RemapUtils.cc

                          

add_amanzi_library(operators 
  SOURCE ${operators_src_files}
  HEADERS ${operators_inc_files}
  LINK_LIBS ${operators_link_libs})


if (BUILD_TESTS)
    # Add UnitTest includes
    include_directories(${UnitTest_INCLUDE_DIRS})

    # Copy test directory files if an out of source build
    if (NOT (${OPERATORS_SOURCE_DIR} EQUAL ${OPERATORS_BINARY_DIR}))
        file(GLOB DataFiles "${OPERATORS_SOURCE_DIR}/test/*.xml"
                            "${OPERATORS_SOURCE_DIR}/test/*.exo")
        file(COPY ${DataFiles} DESTINATION ${OPERATORS_BINARY_DIR}/test/)
    endif()

    message(DEBUG "KokkosKernels_LIBRARIES = ${KokkosKernels_LIBRARIES}")
    set(operators_test_link_libs geometry atk mesh mesh_factory  data_structures operators whetstone output atk ${KokkosKernels_LIBRARIES})
#mesh_logical
    
    # Tests: diffusion
    # add_amanzi_test(operators_laplace_beltrami operators_laplace_beltrami
    #                 KIND unit
    #                 SOURCE test/Main.cc test/operator_laplace_beltrami.cc 
    #                 LINK_LIBS ${operators_test_link_libs})

    # add_amanzi_test(operators_parabolic_surface operators_parabolic_surface
    #                 KIND unit
    #                 SOURCE test/Main.cc test/operator_parabolic_surface.cc
    #                 LINK_LIBS ${operators_test_link_libs})

    # add_amanzi_test(operators_marshak operators_marshak
    #                 KIND unit
    #                 NPROCS 3
    #                 SOURCE test/Main.cc test/operator_marshak.cc
    #                        test/operator_marshak_logical.cc
    #                 LINK_LIBS ${operators_test_link_libs})

    # add_amanzi_test(operators_nonlinear_surface operators_nonlinear_surface
    #                 KIND unit
    #                 SOURCE test/Main.cc test/operator_nonlinear_surface.cc
    #                 LINK_LIBS ${operators_test_link_libs})

    add_amanzi_test(operators_diffusion_discretization operators_diffusion_discretization
                    KIND unit
                    SOURCE test/Main.cc test/operator_diffusion_discretizations.cc
                    LINK_LIBS ${UnitTest_LIBRARIES} ${operators_test_link_libs})

    add_amanzi_test(operators_diffusion_disc_3pe operators_diffusion_discretization NPROCS 3 KIND unit)
                  
    # add_amanzi_test(operators_diffusion operators_diffusion
    #                 KIND unit
    #                 SOURCE test/Main.cc test/operator_diffusion_nodal.cc
    #                                     test/operator_diffusion_divk.cc
    #                                     test/operator_diffusion_nonsymmetric.cc
    #                                     test/operator_diffusion_nlfv.cc
    #                                     test/operator_diffusion_curved.cc
    #                                     test/operator_diffusion_edges.cc test/MFD3D_Diffusion_Edge.cc
    #                 LINK_LIBS ${operators_test_link_libs})

    # add_amanzi_test(operators_diffusion_3pe operators_diffusion NPROCS 3 KIND unit)

    # add_amanzi_test(operators_diffusion_consistent_face operators_diffusion_consistent_face
    #                 KIND unit
    #                 SOURCE test/Main.cc test/operator_diffusion_consistent_face.cc
    #                 LINK_LIBS ${operators_test_link_libs})

    # add_amanzi_test(operators_diffusion_consistent_face_3pe operators_diffusion_consistent_face 
    #                 NPROCS 3 KIND unit)

    # add_amanzi_test(operators_diffusion_dg operators_diffusion_dg
    #                 KIND unit
    #                 SOURCE test/Main.cc test/operator_diffusion_dg.cc
    #                 LINK_LIBS ${operators_test_link_libs})

    # add_amanzi_test(operators_diffusion_fractured_matrix operators_diffusion_fractured_matrix
    #                 KIND unit
    #                 NPROCS 2
    #                 SOURCE test/Main.cc test/operator_diffusion_fractured_matrix.cc
    #                 LINK_LIBS ${operators_test_link_libs})

    # add_amanzi_test(operators_commute operators_commute
    #                 KIND unit
    #                 SOURCE test/Main.cc test/operator_commute.cc
    #                 LINK_LIBS ${operators_test_link_libs})
    
    # # Test: advection-diffusion 
    # add_amanzi_test(operators_advection_diffusion operators_advection_diffusion
    #                 KIND unit
    #                 SOURCE test/Main.cc test/operator_advdiff.cc
    #                 LINK_LIBS ${operators_test_link_libs})

    # # Test: matrices for surface and fractures
    # add_amanzi_test(operators_surface operators_surface
    #                 KIND unit
    #                 SOURCE test/Main.cc test/operator_advdiff_surface.cc
    #                 LINK_LIBS ${operators_test_link_libs})

    # add_amanzi_test(operators_dfn operators_dfn
    #                 KIND unit
    #                 SOURCE test/Main.cc test/operator_advection_dfn.cc
    #                                     test/operator_diffusion_dfn.cc
    #                 LINK_LIBS ${operators_test_link_libs})

    # add_amanzi_test(operators_dfn_4pe operators_dfn NPROCS 4 KIND unit)

    # # Test: electromagnetics 
    # add_amanzi_test(operators_electromagnetics operators_electromagnetics
    #                 KIND unit
    #                 SOURCE test/Main.cc test/AnalyticElectromagneticsBase.cc
    #                                     test/operator_electromagnetics.cc 
    #                                     test/operator_magnetic_diffusion.cc
    #                 LINK_LIBS ${operators_test_link_libs})

    # # Test: utils
    # add_amanzi_test(operators_upwind operators_upwind
    #                 KIND unit
    #                 SOURCE test/Main.cc test/operator_upwind.cc
    #                 LINK_LIBS ${operators_test_link_libs})

    # add_amanzi_test(operators_upwind_2pe operators_upwind NPROCS 2 KIND unit)

    # # stability analysis
    # add_amanzi_test(operators_stability operators_stability
    #                 KIND unit
    #                 SOURCE test/Main.cc test/operator_stability.cc
    #                 LINK_LIBS ${operators_test_link_libs})

    # add_amanzi_test(operators_stability_2pe operators_stability NPROCS 2 KIND unit)

    # # convergence rates
    # add_amanzi_test(operators_diffusion_convergence operators_diffusion_convergence
    #                 KIND unit
    #                 SOURCE test/Main.cc test/operator_diffusion_convergence.cc
    #                 LINK_LIBS ${operators_test_link_libs})
    # add_amanzi_test(operators_coupled_diffusion operators_coupled_diffusion
    #                 KIND unit
    #                 SOURCE test/Main.cc test/operator_coupled_diffusion.cc
    #                 LINK_LIBS ${operators_test_link_libs})
    # add_amanzi_test(operators_coupled_diffusion_2pe operators_coupled_diffusion NPROCS 2 KIND unit)

    
    # # Test: diffusion with gravity
    # add_amanzi_test(operators_gravity operators_gravity
    #                 KIND unit
    #                 SOURCE test/Main.cc test/operator_gravity.cc
    #                 LINK_LIBS ${operators_test_link_libs})

    # # Test: elasticity methods
    # add_amanzi_test(operators_elasticity operators_elasticity
    #                 KIND unit
    #                 SOURCE test/Main.cc test/operator_elasticity.cc
    #                 LINK_LIBS ${operators_test_link_libs})

    # add_amanzi_test(operators_stokes operators_stokes
    #                 KIND unit
    #                 SOURCE test/Main.cc test/operator_stokes.cc
    #                 LINK_LIBS ${operators_test_link_libs})

    # # reconstruction tests
    # add_amanzi_test(operators_reconstruction_limiters operators_reconstruction_limiters
    #                 KIND unit
    #                 SOURCE test/Main.cc test/operator_reconstruction.cc
    #                                     test/operator_limiters.cc
    #                                     test/operator_limiters_dg.cc
    #                 LINK_LIBS ${operators_test_link_libs})

    # add_amanzi_test(operators_boundary_flux_solver operators_boundary_flux_solver
    #                 KIND unit
    #                 SOURCE test/Main.cc test/operator_boundary_flux_solver.cc
    #                 LINK_LIBS ${operators_test_link_libs})

    # add_amanzi_test(operators_reconstruction_limiters_2pe operators_reconstruction_limiters
    #                 NPROCS 2
    #                 KIND unit)

    # # remap tests
    # add_amanzi_test(operators_remap operators_remap
    #                 KIND unit
    #                 NPROCS 2
    #                 SOURCE test/Main.cc test/operator_remap.cc
    #                 LINK_LIBS ${operators_test_link_libs})

    # add_amanzi_test(operators_remap_curved operators_remap_curved
    #                 KIND unit
    #                 SOURCE test/Main.cc test/operator_remap_curved.cc
    #                 LINK_LIBS ${operators_test_link_libs})

    # # high-order diffusion test
    # add_amanzi_test(operators_diffusion_high_order operators_diffusion_high_order
    #                 KIND unit
    #                 SOURCE test/Main.cc test/operator_diffusion_high_order.cc
    #                 LINK_LIBS ${operators_test_link_libs})

    # # DG advection test
    # add_amanzi_test(operators_advection_dg operators_advection_dg
    #                 KIND unit
    #                 NPROCS 3
    #                 SOURCE test/Main.cc test/operator_advection_dg.cc
    #                                     test/operator_advection_dg_transient
    #                 LINK_LIBS ${operators_test_link_libs})

    # # Test: mini PDEs
    # add_amanzi_test(operators_mini1D operators_mini1D
    #                 KIND unit
    #                 SOURCE test/Main.cc test/operator_mini1D.cc
    #                                     test/operator_mini1D_nonlinear.cc
    #                 LINK_LIBS ${operators_test_link_libs})

    # # DG advection test
    # add_amanzi_test(operators_advection_dg operators_advection_dg
    #                 KIND unit
    #                 NPROCS 3
    #                 SOURCE test/Main.cc test/operator_advection_dg.cc
    #                                     test/operator_advection_dg_transient
    #                 LINK_LIBS ${operators_test_link_libs})

endif()


