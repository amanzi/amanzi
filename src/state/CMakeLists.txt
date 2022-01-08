# -*- mode: cmake -*-

#
#  Amanzi
#    State class
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)
include(LibraryManager)

#
# Define a project name
# After this command the following varaibles are defined
#   STATE_SOURCE_DIR
#   STATE_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(STATE)

# Amanzi include directories
include_directories(${ATK_SOURCE_DIR})
include_directories(${DATA_STRUCTURES_SOURCE_DIR})
include_directories(${DBC_SOURCE_DIR})
include_directories(${DBG_SOURCE_DIR})
include_directories(${FUNCS_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${MESH_MSTK_SOURCE_DIR})
include_directories(${MFUNCS_SOURCE_DIR})
include_directories(${OPERATORS_SOURCE_DIR})
include_directories(${OUTPUT_SOURCE_DIR})
include_directories(${SOLVERS_SOURCE_DIR})
include_directories(${STATE_SOURCE_DIR}/data)
include_directories(${WHETSTONE_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${ASCEMIO_INCLUDE_DIR})
include_directories(${SEACAS_DIR}/include)
if(ENABLE_Silo)
  include_directories(${Silo_INCLUDE_DIRS})
  add_definitions(-DENABLE_Silo)
endif(ENABLE_Silo)


# This state library will move to a new location 
# once Markus has created the data manager.
# I use global properties since all the PROJECT_NAME
# variables in this directory are tied to the MPC.
# These properties will be replaced with the appropriate
# directory project name once it is moved.
#
# Library: state
#
set(state_inc_files 
  Checkpoint.hh
  chemistry_data.hh 
  data/Data.hh
  data/Data_Helpers.hh
  data/Data_Impl.hh
  data/DataFactory.hh
  data/DataFactory_Impl.hh
  data/Record.hh
  Evaluator.hh
  Evaluator_Factory.hh
  Evaluator_OperatorApply.hh
  Evaluator_PDE_Diffusion.hh
  EvaluatorCellVolume.hh
  EvaluatorDeformingCellVolume.hh
  EvaluatorIndependent.hh
  EvaluatorIndependentFromFile.hh
  EvaluatorIndependentFunction.hh
  EvaluatorMultiplicativeReciprocal.hh
  EvaluatorPrimary.hh
  EvaluatorSecondaryMonotype.hh
  EvaluatorSecondaryMonotypeFromFunction.hh
  IO.hh
  Observable.hh
  StateDefs.hh
  State.hh
  Visualization.hh
  VisualizationDomainSet.hh
  UnstructuredObservations.hh
  Debugger.hh)

add_amanzi_library(state 
                   SOURCE State.cc
                          Checkpoint.cc
                          chemistry_data.cc
                          Evaluator.cc
                          Evaluator_Factory.cc
                          Evaluator_PDE_Diffusion.cc
                          Evaluator_OperatorApply.cc
                          EvaluatorCellVolume.cc
                          EvaluatorDeformingCellVolume.cc
                          EvaluatorIndependent.cc
                          EvaluatorIndependentFromFile.cc
                          EvaluatorIndependentFunction.cc
                          EvaluatorMultiplicativeReciprocal.cc
                          EvaluatorPrimary.cc
                          EvaluatorSecondary.cc
                          EvaluatorSecondaryMonotype.cc
                          EvaluatorSecondaryMonotypeFromFunction.cc
                          IO.cc
                          Observable.cc
                          UnstructuredObservations.cc
                          Visualization.cc
                          VisualizationDomainSet.cc
                          data/Record.cc
                          data/RecordSet.cc
                          data/Data_Helpers.cc
                          Debugger.cc
                   HEADERS ${state_inc_files}
		   LINK_LIBS operators data_structures functions mesh_functions output mesh atk
		   ${Boost_LIBRARIES} ${SEACAS_LIBRARIES})
set_property(GLOBAL PROPERTY STATE_SOURCE_DIR ${MPC_SOURCE_DIR})
set_property(GLOBAL PROPERTY STATE_BINARY_DIR ${MPC_BINARY_DIR})
if (APPLE AND BUILD_SHARED_LIBS)
  set_target_properties(state PROPERTIES LINK_FLAGS "-Wl,-undefined,dynamic_lookup")
endif()			  

#
# Install Target
#
add_install_include_file(${state_inc_files})

if (BUILD_TESTS)
    # Add UnitTest includes
    include_directories(${UnitTest_INCLUDE_DIRS})

    # Add State include directory. 
    include_directories(${STATE_SOURCE_DIR})

    # Mesh factory includes
    include_directories(${MESH_FACTORY_SOURCE_DIR})

    # Copy test subdirectory for out of source builds
    if (NOT ("${STATE_SOURCE_DIR}" STREQUAL "${STATE_BINARY_DIR}"))
      file(GLOB DataFiles "${STATE_SOURCE_DIR}/test/*.xml"
                          "${STATE_SOURCE_DIR}/test/*.exo"
                          "${STATE_SOURCE_DIR}/test/*.gold")
      file(COPY ${DataFiles} DESTINATION ${STATE_BINARY_DIR}/test/)
    endif()
    
    add_amanzi_test(state_restart state_restart
                    KIND int
                    SOURCE
	  	          test/Main.cc
	 		  test/state_restart.cc
	 		  test/state_vis.cc
                    LINK_LIBS state mesh_factory mesh data_structures operators mesh_functions functions 
                              ${UnitTest_LIBRARIES} ${NOX_LIBRARIES} ${HDF5_LIBRARIES})

    add_amanzi_test(state_restart2 state_restart2
                    KIND int
                    SOURCE
	  	          test/Main.cc
	 		  test/state_restart2.cc
                    NPROCS 2
                    LINK_LIBS state mesh_factory mesh data_structures operators mesh_functions functions 
                              ${UnitTest_LIBRARIES} ${NOX_LIBRARIES} ${HDF5_LIBRARIES})
                            
    add_amanzi_test(state_copy state_copy
                    KIND int
                    SOURCE test/Main.cc test/state_copy.cc
                    LINK_LIBS state mesh_factory mesh data_structures operators mesh_functions functions 
                              ${UnitTest_LIBRARIES} ${NOX_LIBRARIES} ${HDF5_LIBRARIES})

    add_amanzi_test(state_dag state_dag
                    KIND int
                    SOURCE test/Main.cc test/state_dag.cc
                    LINK_LIBS state mesh_factory mesh data_structures operators mesh_functions functions 
                              ${UnitTest_LIBRARIES} ${NOX_LIBRARIES} ${HDF5_LIBRARIES})

    add_amanzi_test(state_functional_evaluators state_functional_evaluators
                    KIND int
                    SOURCE test/Main.cc
	 		   test/state_functional_evaluators.cc
                    LINK_LIBS state mesh_factory mesh data_structures operators mesh_functions functions 
                              ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})

    add_amanzi_test(state_init state_init
                    KIND int
                    SOURCE test/Main.cc
	 		   test/state_init.cc
                    LINK_LIBS state mesh_factory mesh data_structures operators mesh_functions functions 
                              ${UnitTest_LIBRARIES} ${NOX_LIBRARIES} ${HDF5_LIBRARIES})

    add_amanzi_test(state_observations state_observations
                    KIND int
                    SOURCE test/Main.cc
	 		   test/state_observations.cc
                    LINK_LIBS state mesh_factory mesh data_structures operators mesh_functions functions 
                              ${UnitTest_LIBRARIES} ${NOX_LIBRARIES} ${HDF5_LIBRARIES})
    add_amanzi_test(state_observations_np2 state_observations NPROCS 2 KIND unit)

    add_amanzi_test(state_data_factory state_data_factory
                    KIND unit
                    SOURCE test/Main.cc 
                           test/state_data_factory.cc
                    LINK_LIBS error_handling atk mesh operators data_structures mesh_functions state
                              ${state_tpl_list} ${UnitTest_LIBRARIES})

    add_amanzi_test(state_record state_record
                    KIND unit
                    SOURCE test/Main.cc
                           test/state_record.cc
                    LINK_LIBS atk mesh operators data_structures mesh_functions state
                              ${state_tpl_list} ${UnitTest_LIBRARIES})

    add_amanzi_test(state_mesh state_mesh
                    KIND unit
                    SOURCE test/Main.cc
                           test/state_mesh.cc
                    LINK_LIBS atk mesh_factory mesh operators data_structures mesh_functions state
                              ${state_tpl_list} ${UnitTest_LIBRARIES})

    add_amanzi_test(state_vis state_vis
                    KIND int
                    SOURCE test/Main.cc
                           test/state_vis.cc
                    LINK_LIBS state mesh_factory mesh operators data_structures mesh_functions functions
                              ${state_tpl_list} ${UnitTest_LIBRARIES})

    add_amanzi_test(state_container state_container
                    KIND unit
                    SOURCE test/Main.cc
                           test/state_container.cc
                    LINK_LIBS atk mesh_factory mesh operators data_structures mesh_functions operators state
                              ${state_tpl_list} ${UnitTest_LIBRARIES})

    add_amanzi_test(state_extensibility state_extensibility
                    KIND int
                    SOURCE test/Main.cc
                           test/state_extensibility.cc
                    LINK_LIBS state mesh_factory mesh operators data_structures mesh_functions functions 
                              ${state_tpl_list} ${UnitTest_LIBRARIES})
    
    add_amanzi_test(state_evaluators state_evaluators
                    KIND int
                    SOURCE test/Main.cc
                           test/state_evaluators.cc
                    LINK_LIBS state mesh_factory mesh operators data_structures mesh_functions functions 
                              ${UnitTest_LIBRARIES})

    add_amanzi_test(state_evaluators_cv state_evaluators_cv
                    KIND int
                    SOURCE test/Main.cc
                           test/state_evaluators_cv.cc
                    LINK_LIBS state mesh_factory mesh operators data_structures mesh_functions functions 
                              ${UnitTest_LIBRARIES})

    add_amanzi_test(state_evaluators_op state_evaluators_op
                    KIND int
                    SOURCE test/Main.cc
                           test/state_evaluators_op.cc
                    LINK_LIBS state mesh_factory mesh operators data_structures mesh_functions functions 
                              ${UnitTest_LIBRARIES})

    add_amanzi_test(state_evaluators_op_inverse state_evaluators_op_inverse
                    KIND int
                    SOURCE test/Main.cc
                           test/state_evaluators_op_inverse.cc
                    LINK_LIBS state mesh_factory mesh operators data_structures mesh_functions functions 
                              ${UnitTest_LIBRARIES})

endif()

