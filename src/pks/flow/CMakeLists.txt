#
#  Amanzi
#   Flow process kernel
#
# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)

# Define a project name
# After this command the following varaibles are defined
#   FLOW_SOURCE_DIR
#   FLOW_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(FLOW)

# External (TPL) include directories
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${Teuchos_INCLUDE_DIRS})

# Amanzi include directories
include_directories(${ATK_SOURCE_DIR})
include_directories(${DATA_STRUCTURES_SOURCE_DIR})
include_directories(${EOS_SOURCE_DIR})
include_directories(${EVALUATORS_SOURCE_DIR})
include_directories(${FUNCS_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${MFUNCS_SOURCE_DIR})
include_directories(${MESH_MSTK_SOURCE_DIR})
include_directories(${OUTPUT_SOURCE_DIR})
include_directories(${OPERATORS_SOURCE_DIR})
include_directories(${PKS_SOURCE_DIR})
include_directories(${SOLVERS_SOURCE_DIR})
include_directories(${STATE_SOURCE_DIR})
include_directories(${TIME_INTEGRATION_SOURCE_DIR})
include_directories(${WHETSTONE_SOURCE_DIR})

#
# Library: flow
#
set(flow_src_files Flow_PK.cc 
                   Flow_SeepageFace.cc FlowBoundaryFunction.cc
                   Flow_VandV.cc  
                   Darcy_PK.cc Darcy_TI.cc Darcy_LinearOperator.cc
                   Richards_PK.cc Richards_TI.cc Richards_WRM.cc
                   Richards_Picard.cc Richards_LinearOperator.cc Richards_CNLS.cc
                   ApertureDarcyEvaluator.cc
                   PermeabilityModelPartition.cc PermeabilityEvaluator.cc
                   RelPermEvaluator.cc
                   WRM_vanGenuchten.cc WRM_BrooksCorey.cc WRM_linear.cc WRM_tabular.cc
                   WRMPartition.cc WRMFactory.cc WRMEvaluator.cc 
                   WaterStorage.cc WaterStorageDarcy.cc
                   DarcyVelocityEvaluator.cc
                   VolumetricFlowRateEvaluator.cc
                   MultiscaleFlowPorosityFactory.cc MultiscaleFlowPorosityPartition.cc
                   MultiscaleFlowPorosity_DPM.cc
                   MultiscaleFlowPorosity_GDPM.cc
                   FracturePermModelPartition.cc FracturePermModelEvaluator.cc
                   )
set(flow_tpl_libs ${Teuchos_LIBRARIES} ${Epetra_LIBRARIES} ${SEACAS_LIBRARIES})
if (Ifpack_ENABLE_HYPRE)
  list(APPEND flow_tpl_libs ${Ifpack_LIBRARIES} ${HYPRE_LIBRARY})
endif()
add_amanzi_library(flow SOURCE ${flow_src_files}
                   LINK_LIBS pks eos evaluators operators state functions
                             mesh_functions time_integration data_structures
                             ${flow_tpl_libs})

#
# Install Targets
#
set(flow_inc_files 
  Darcy_PK.hh
  FlowDefs.hh
  FlowBoundaryFunction.hh
  FlowSourceFunction.hh
  Flow_PK.hh
  FlowTypeDefs.hh
  FracturePermModel.hh
  FracturePermModel_Constant.hh
  FracturePermModel_CubicLaw.hh
  FracturePermModel_Linear.hh
  FracturePermModelPartition.hh
  FracturePermModelEvaluator.hh
  models_flow_reg.hh
  MultiscaleFlowPorosity.hh
  MultiscaleFlowPorosity_DPM.hh
  MultiscaleFlowPorosity_GDPM.hh
  MultiscaleFlowPorosityFactory.hh
  MultiscaleFlowPorosityPartition.hh
  Richards_PK.hh
  VolumetricFlowRateEvaluator.hh
  WRM.hh
  WRM_BrooksCorey.hh
  WRM_linear.hh
  WRM_saturated.hh
  WRM_tabular.hh
  WRM_vanGenuchten.hh
  WRMEvaluator.hh
  WRMFactory.hh
  WRMPartition.hh)

add_install_include_file(${flow_inc_files})


if (BUILD_TESTS) 
    # Add UnitTest include directory
    include_directories(${UnitTest_INCLUDE_DIRS})
    include_directories(${MESH_FACTORY_SOURCE_DIR})

    # Copy test directory files if any out of source build
    if (NOT (${FLOW_SOURCE_DIR} EQUAL ${FLOW_BINARY_DIR}) )
        # execute_process(COMMAND ${CMAKE_COMMAND} -E 
        #   copy_directory ${FLOW_SOURCE_DIR}/test ${FLOW_BINARY_DIR}/test) 
        file(GLOB DataFiles "${FLOW_SOURCE_DIR}/test/*.xml"
                            "${FLOW_SOURCE_DIR}/test/*.exo"
                            "${FLOW_SOURCE_DIR}/test/*.h5")
        file(COPY ${DataFiles} DESTINATION ${FLOW_BINARY_DIR}/test/)
    endif()

    # Add the flow directory to the include paths
    include_directories(${FLOW_SOURCE_DIR})

    set(amanzi_libs geometry mesh mesh_factory state flow pks operators whetstone data_structures)

    # Test: transient Darcy flow 
    add_amanzi_test(flow_darcy_transient flow_darcy_transient
                    KIND int
                    SOURCE test/Main.cc test/flow_darcy_transient.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})
    add_amanzi_test(flow_darcy_transient_parallel flow_darcy_transient NPROCS 2 KIND unit)
 
    # Test: transient Darcy flow 
    add_amanzi_test(flow_darcy_dual2D flow_darcy_dual2D
                    KIND int
                    SOURCE test/Main.cc test/flow_darcy_dual2D.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})
 
    # Test: transient Darcy flow with wells
    add_amanzi_test(flow_darcy_well flow_darcy_well
                    KIND int
                    SOURCE test/Main.cc test/flow_darcy_well.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})

    # Test: transient Darcy flow with sources
    add_amanzi_test(flow_darcy_source flow_darcy_source
                    KIND int
                    SOURCE test/Main.cc test/flow_darcy_source.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})

    # Test: convergence analysis
    add_amanzi_test(flow_richards_convergence flow_richards_convergence
                    KIND int
                    SOURCE test/Main.cc test/flow_richards_convergence.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})

    # Test: convergence analysis on random meshes
    add_amanzi_test(flow_richards_random flow_richards_random
                    KIND int
                    SOURCE test/Main.cc test/flow_richards_random.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})

    # Test: 2D Richards
    add_amanzi_test(flow_richards_2D flow_richards_2D
                    KIND unit
                    SOURCE test/Main.cc test/flow_richards_2D.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})

    add_amanzi_test(flow_richards_seepage flow_richards_seepage
                    KIND unit
                    SOURCE test/Main.cc test/flow_richards_seepage.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})

    add_amanzi_test(flow_richards_porosity flow_richards_porosity
                    KIND unit
                    SOURCE test/Main.cc test/flow_richards_porosity.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})

    # Test: 2D Richards TPFA discretization test seepage BC
    add_amanzi_test(flow_richards_seepage_tpfa flow_richards_seepage_tpfa
                    KIND unit
                    SOURCE test/Main.cc test/flow_richards_seepage_tpfa.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})

    # Test: pseudo-1D Richards
    add_amanzi_test(flow_richards_bc_cribs flow_richards_bc_cribs
                    KIND unit
                    SOURCE test/Main.cc test/flow_richards_bc_cribs.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})

    # Test: Newton
    add_amanzi_test(flow_richards_newton_tpfa flow_richards_newton_tpfa
                    KIND unit
                    SOURCE test/Main.cc test/flow_richards_newton_tpfa.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})
 
    # Test: various sources
    add_amanzi_test(flow_richards_source flow_richards_source
                    KIND int
                    SOURCE test/Main.cc test/flow_richards_source.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})

    # Test: miscaleneous Darcy flow routines
    add_amanzi_test(flow_darcy_misc flow_darcy_misc 
                    KIND int
                    SOURCE test/Main.cc test/flow_darcy_misc.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})

    add_amanzi_test(flow_darcy_misc_2pe flow_darcy_misc NPROCS 2 KIND unit)
    add_amanzi_test(flow_darcy_misc_4pe flow_darcy_misc NPROCS 4 KIND unit)

    # Test: model tests
    add_amanzi_test(flow_wrm flow_wrm
                    KIND int
                    SOURCE test/Main.cc test/flow_vanGenuchten.cc
                                        test/flow_wrm_tabular.cc
                                        test/flow_wrm_solvers.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})

    # Test: flow-bc-factory
    add_amanzi_test(flow_bc_factory flow_bc_factory
                    KIND unit
                    SOURCE test/Main.cc test/flow_bc_factory.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})

    # Test: flow-boundary solver test
    add_amanzi_test(flow_richards_boundary_solver flow_richards_boundary_solver
                    KIND unit
                    SOURCE test/Main.cc test/flow_richards_boundary_solver.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})
 
    # Test: 3D Richards
    set(amanzi_libs geometry mesh mesh_factory flow pks operators state whetstone solvers output data_structures)
    
    add_amanzi_test(flow_richards_tensor flow_richards_tensor
                    KIND unit
                    SOURCE test/Main.cc test/flow_richards_tensor.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})

    add_amanzi_test(flow_richards_3D flow_richards_3D
                    KIND unit
                    SOURCE test/Main.cc test/flow_richards_3D.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})

    # Test: Multiscale
    add_amanzi_test(flow_multiscale flow_multiscale
                    KIND unit
                    SOURCE test/Main.cc test/flow_multiscale.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})

    if(ENABLE_MESH_MSTK)
      # Test: Fractures
      add_amanzi_test(flow_darcy_fractures flow_darcy_fractures 
                      KIND int
                      SOURCE test/Main.cc test/flow_darcy_fractures.cc
                      LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES}
                      NPROCS 4)

      add_amanzi_test(flow_darcy_fractures_source flow_darcy_fractures_source
                      KIND int
                      SOURCE test/Main.cc test/flow_darcy_fractures_source.cc
                      LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})

      add_amanzi_test(flow_richards_fractures_saturated flow_richards_fractures_saturated 
                      KIND int
                      SOURCE test/Main.cc test/flow_richards_fractures_saturated.cc
                      LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})

      add_amanzi_test(flow_richards_fractures flow_richards_fractures
                      KIND int
                      SOURCE test/Main.cc test/flow_richards_fractures.cc
                      LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})
    endif()
                  
endif()

