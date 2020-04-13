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
include_directories(${OPERATORS_SOURCE_DIR})
include_directories(${WHETSTONE_SOURCE_DIR})
include_directories(${DATA_STRUCTURES_SOURCE_DIR})
include_directories(${DBC_SOURCE_DIR})
include_directories(${DBG_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${OUTPUT_SOURCE_DIR})
include_directories(${FUNCS_SOURCE_DIR})
include_directories(${MFUNCS_SOURCE_DIR})
include_directories(${SOLVERS_SOURCE_DIR})
include_directories(${STATE_SOURCE_DIR}/data)
include_directories(${STATE_SOURCE_DIR}/io)

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${ASCEMIO_INCLUDE_DIR})
include_directories(${ExodusII_INCLUDE_DIR})
if(ENABLE_Silo)
  include_directories(${Silo_INCLUDE_DIRS})
  add_definitions(-DENABLE_Silo)
endif(ENABLE_Silo)

#
# Library: state
#
set(state_inc_files 
  data/Data.hh
  data/Data_Initializers.hh
  data/Data_Impl.hh
  data/Record.hh
  data/RecordSet.hh
  data/DataFactory.hh
  data/DataFactory_Impl.hh
  evaluator/Evaluator.hh
  evaluator/Evaluator_Factory.hh
  evaluator/EvaluatorIndependent.hh
  evaluator/EvaluatorIndependentFromFile.hh
  evaluator/EvaluatorIndependentFunction.hh
  evaluator/EvaluatorIndependentConstant.hh
  evaluator/EvaluatorIndependentTensorFunction.hh
  evaluator/EvaluatorPrimary.hh
  evaluator/EvaluatorSecondary.hh
  evaluator/EvaluatorSecondaryMonotype.hh
  evaluator/EvaluatorSecondaryMonotypeFromFunction.hh
  evaluator/Evaluator_OperatorApply.hh
  evaluator/Evaluator_PDE_Diffusion.hh
  evaluator/Evaluator_PDE_Accumulation.hh
  evaluator/Evaluator_BCs.hh
  evaluator/Evaluator_BCsFunction.hh
  evaluator/EvaluatorSecondaryMonotypeAdditive.hh
  evaluator/EvaluatorSecondaryMonotypeMultiplicative.hh
  evaluator/EvaluatorCellVolume.hh
#  evaluator/cell_volume_evaluator.hh
#  evaluator/deforming_cell_volume_evaluator.hh
#  evaluator/rank_evaluator.hh
  io/Checkpoint.hh
  io/Debugger.hh
  io/Observable.hh
  io/UnstructuredObservations.hh
  io/Visualization.hh
  StateDefs.hh
  State.hh
)


set(state_src_files 
  data/Data_Initializers.cc
  data/Record.cc
  data/RecordSet.cc
  evaluator/Evaluator.cc
  evaluator/Evaluator_Factory.cc
  evaluator/EvaluatorIndependent.cc
  evaluator/EvaluatorIndependentFromFile.cc
  evaluator/EvaluatorIndependentFunction.cc
  evaluator/EvaluatorIndependentConstant.cc
  evaluator/EvaluatorIndependentTensorFunction.cc
  evaluator/EvaluatorPrimary.cc
  evaluator/EvaluatorSecondary.cc
  evaluator/EvaluatorSecondaryMonotype.cc
  evaluator/EvaluatorSecondaryMonotypeFromFunction.cc
  evaluator/Evaluator_OperatorApply.cc
  evaluator/Evaluator_PDE_Diffusion.cc
  evaluator/Evaluator_PDE_Accumulation.cc
  evaluator/Evaluator_BCs.cc
  evaluator/EvaluatorCellVolume.cc
#  evaluator/cell_volume_evaluator.cc
#  evaluator/deforming_cell_volume_evaluator.cc
#  evaluator/rank_evaluator.cc
  io/Debugger.cc
  io/Checkpoint.cc
  io/Visualization.cc
  io/Observable.cc
  io/UnstructuredObservations.cc
  State.cc
  )

set(state_link_libs
  error_handling
  atk
  mesh
  data_structures
  mesh_functions
  output
  operators
  ${Boost_LIBRARIES}
  )

#
# Library: state
#
add_amanzi_library(state 
                   SOURCE ${state_src_files}
                   HEADERS ${state_inc_files}
                   LINK_LIBS ${state_link_libs}
                   )

#
# Install header files
#
add_install_include_file(${state_inc_files})

#
# Register with factories
#
register_evaluator_with_factory(HEADERFILE evaluator/EvaluatorPrimary_reg.hh LISTNAME REGISTER_AMANZI_STATE_EVALUATORS)
register_evaluator_with_factory(HEADERFILE evaluator/EvaluatorIndependentFunction_reg.hh LISTNAME REGISTER_AMANZI_STATE_EVALUATORS)
register_evaluator_with_factory(HEADERFILE evaluator/EvaluatorIndependentTensorFunction_reg.hh LISTNAME REGISTER_AMANZI_STATE_EVALUATORS)
register_evaluator_with_factory(HEADERFILE evaluator/EvaluatorIndependentFromFile_reg.hh LISTNAME REGISTER_AMANZI_STATE_EVALUATORS)
register_evaluator_with_factory(HEADERFILE evaluator/EvaluatorIndependentConstant_reg.hh LISTNAME REGISTER_AMANZI_STATE_EVALUATORS)
register_evaluator_with_factory(HEADERFILE evaluator/Evaluator_OperatorApply_reg.hh LISTNAME REGISTER_AMANZI_STATE_EVALUATORS)
register_evaluator_with_factory(HEADERFILE evaluator/Evaluator_PDE_Diffusion_reg.hh LISTNAME REGISTER_AMANZI_STATE_EVALUATORS)
register_evaluator_with_factory(HEADERFILE evaluator/Evaluator_PDE_Accumulation_reg.hh LISTNAME REGISTER_AMANZI_STATE_EVALUATORS)
register_evaluator_with_factory(HEADERFILE evaluator/Evaluator_BCs_reg.hh LISTNAME REGISTER_AMANZI_STATE_EVALUATORS)
register_evaluator_with_factory(HEADERFILE evaluator/Evaluator_BCsFunction_reg.hh LISTNAME REGISTER_AMANZI_STATE_EVALUATORS)
register_evaluator_with_factory(HEADERFILE evaluator/EvaluatorSecondaryMonotypeAdditive_reg.hh LISTNAME REGISTER_AMANZI_STATE_EVALUATORS)
register_evaluator_with_factory(HEADERFILE evaluator/EvaluatorSecondaryMonotypeMultiplicative_reg.hh LISTNAME REGISTER_AMANZI_STATE_EVALUATORS)
register_evaluator_with_factory(HEADERFILE evaluator/EvaluatorSecondaryMonotypeFromFunction_reg.hh LISTNAME REGISTER_AMANZI_STATE_EVALUATORS)
register_evaluator_with_factory(HEADERFILE evaluator/EvaluatorCellVolume_reg.hh LISTNAME REGISTER_AMANZI_STATE_EVALUATORS)
# register_evaluator_with_factory(HEADERFILE evaluator/rank_evaluator_reg.hh          LISTNAME REGISTER_AMANZI_STATE_EVALUATORS)
# register_evaluator_with_factory(HEADERFILE evaluator/deforming_cell_volume_evaluator_reg.hh LISTNAME REGISTER_AMANZI_STATE_EVALUATORS)

if (REGENERATE_REGISTRATIONS)
  generate_evaluators_registration_header(HEADERFILE state_evaluators_registration.hh LISTNAME REGISTER_AMANZI_STATE_EVALUATORS  INSTALL True)
endif()

#
# Build tests
#
if (BUILD_TESTS)
    # Add UnitTest includes
    include_directories(${UnitTest_INCLUDE_DIRS})

    # other libraries
    include_directories(${STATE_SOURCE_DIR})
    include_directories(${MESH_FACTORY_SOURCE_DIR})

    set(state_test_link_libs state mesh_factory ${UnitTest_LIBRARIES})
    
    # Copy test subdirectory for out of source builds
    if (NOT ("${STATE_SOURCE_DIR}" STREQUAL "${STATE_BINARY_DIR}"))
        execute_process(COMMAND ${CMAKE_COMMAND} -E 
          copy_directory ${STATE_SOURCE_DIR}/test ${STATE_BINARY_DIR}/test) 
    endif()

    # Tests the DataFactory object, a part of state's container capabilities.
    add_amanzi_test(state_data_factory state_data_factory
                    KIND unit
                    SOURCE test/Main.cc test/state_data_factory.cc
                    LINK_LIBS ${state_test_link_libs})

    # Tests the Record and RecordSet objects, meta-data and a part of
    # state's container capabilities.
    add_amanzi_test(state_record state_record
                    KIND unit
                    SOURCE test/Main.cc test/state_record.cc
                    LINK_LIBS ${state_test_link_libs})

    # Tests the ability of state to hold meshes.
    add_amanzi_test(state_mesh state_mesh
                    KIND unit
                    SOURCE test/Main.cc test/state_mesh.cc
                    LINK_LIBS ${state_test_link_libs})

    # Tests the ability of state to hold known data types.
    add_amanzi_test(state_container state_container
                    KIND unit
                    SOURCE test/Main.cc test/state_container.cc
                    LINK_LIBS ${state_test_link_libs})

    # Tests the ability of state to hold user-provided data types.
    add_amanzi_test(state_extensibility state_extensibility
                    KIND int
                    SOURCE test/Main.cc test/state_extensibility.cc
                    LINK_LIBS ${state_test_link_libs})

    # Tests the ability to visualize state.                    
    add_amanzi_test(state_vis state_vis
                    KIND int
                    SOURCE test/Main.cc test/state_vis.cc
                    LINK_LIBS ${state_test_link_libs})

    # Tests checkpoint/restart
    add_amanzi_test(state_restart state_restart
                    KIND int
                    SOURCE test/Main.cc test/state_restart.cc
                    LINK_LIBS ${state_test_link_libs})
    add_amanzi_test(state_restart_np2 state_restart NPROCS 2 KIND unit)

    # Test: TODO: Write test for initializing data from XML.
    
    # Tests the logic of base class evaluators, lazy evaluation, etc.
    add_amanzi_test(state_evaluators_basic state_evaluators_basic
                    KIND int
                    SOURCE test/Main.cc test/state_evaluators_basic.cc
                    LINK_LIBS ${state_test_link_libs})

    # Tests evaluators that are fully resolved without physics,
    # e.g. additive, independent evaluators from functions, etc.
    add_amanzi_test(state_evaluators_functions state_evaluators_functions
                    KIND int
                    SOURCE test/Main.cc test/state_evaluators_functions.cc
                    LINK_LIBS ${state_test_link_libs})

    # Tests EvaluatorSecondaryMonotype's EnsureConsistency, which
    # passes dependencies up and down the dag.
    add_amanzi_test(state_evaluators_secondary_ensureconsistency state_evaluators_secondary_ensureconsistency
                    KIND int
                    SOURCE test/Main.cc test/state_evaluators_secondary_ensureconsistency.cc
                    LINK_LIBS ${state_test_link_libs})

    # Tests a full dag using doubles              
    add_amanzi_test(state_evaluators_dag_double state_evaluators_dag_double
                    KIND int
                    SOURCE test/Main.cc test/state_evaluators_dag_double.cc
                    LINK_LIBS ${state_test_link_libs})

    # Same test but with CompositeVectors and hard-coded evaluators
    add_amanzi_test(state_evaluators_dag_tpetra state_evaluators_dag_tpetra 
          KIND int
          SOURCE test/Main.cc test/state_evaluators_dag_tpetra.cc
          LINK_LIBS ${state_test_link_libs})

    # Same test but with Model and a generic evaluator class
    add_amanzi_test(state_evaluators_dag_models state_evaluators_dag_models
          KIND int
          SOURCE test/Main.cc test/state_evaluators_dag_models.cc
          LINK_LIBS ${state_test_link_libs})

    # Tests a simple physics-based dag using factory -- this is quite
    # close to how users will use State
    #
    # NOTE: this uses a special main to include the static
    # registrations of the test models.
    add_amanzi_test(state_evaluators_factory state_evaluators_factory
          KIND int
          SOURCE test/Main2.cc test/state_evaluators_factory.cc
          LINK_LIBS ${state_test_link_libs})
        
    #-add_amanzi_test(state_evaluators_op state_evaluators_op
    #  KIND int
    #  SOURCE test/Main.cc test/state_evaluators_op.cc
    #  LINK_LIBS ${state_test_link_libs})

    #-add_amanzi_test(state_evaluators_op_inverse state_evaluators_op_inverse
    #  KIND int
    #  SOURCE test/Main.cc test/state_evaluators_op_inverse.cc
    #  LINK_LIBS ${state_test_link_libs})


endif()

