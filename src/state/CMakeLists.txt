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
  data/Data.hh
  data/Data_Helpers.hh
  data/Data_Impl.hh
  data/Record.hh
  data/RecordSet.hh
  data/DataFactory.hh
  data/DataFactory_Impl.hh
  evaluator/Evaluator.hh
  evaluator/Evaluator_Factory.hh
  evaluator/EvaluatorIndependent.hh
  evaluator/EvaluatorIndependentFromFile.hh
  evaluator/EvaluatorIndependentFromFile_reg.hh
  evaluator/EvaluatorIndependentFunction.hh
  evaluator/EvaluatorIndependentFunction_reg.hh
  evaluator/EvaluatorIndependentTensorFunction.hh
  evaluator/EvaluatorIndependentTensorFunction_reg.hh
  evaluator/EvaluatorPrimary.hh
  evaluator/EvaluatorPrimary_reg.hh
  evaluator/EvaluatorSecondary.hh
  evaluator/EvaluatorSecondaryMonotype.hh
  evaluator/EvaluatorSecondaryMonotypeFromFunction.hh
  evaluator/Evaluator_OperatorApply.hh
#  evaluator/Evaluator_PDE_Diffusion.hh
#  evaluator/Evaluator_PDE_Accumulation.hh
  evaluator/Evaluator_BCs.hh
  evaluator/EvaluatorSecondaryMonotypeAdditive.hh
  evaluator/EvaluatorSecondaryMonotypeMultiplicative.hh
  evaluator/EvaluatorCellVolume.hh
#  evaluator/cell_volume_evaluator.hh
#  evaluator/cell_volume_evaluator_reg.hh
#  evaluator/deforming_cell_volume_evaluator.hh
#  evaluator/deforming_cell_volume_evaluator_reg.hh
#  evaluator/rank_evaluator.hh
#  evaluator/rank_evaluator_reg.hh
  io/Checkpoint.hh
  io/Debugger.hh
  io/Observable.hh
  io/UnstructuredObservations.hh
  io/Visualization.hh
  state_evaluators_registration.hh
  StateDefs.hh
  State.hh
)


set(state_tpl_list ${Boost_LIBRARIES} ${Teuchos_LIBRARIES} ${Epetra_LIBRARIES} ${ExodusII_LIBRARIES})
add_amanzi_library(state 
                   SOURCE
                      data/Data_Helpers.cc
		      data/Record.cc
		      data/RecordSet.cc
		      evaluator/Evaluator.cc
		      evaluator/Evaluator_Factory.cc
		      evaluator/EvaluatorIndependent.cc
		      evaluator/EvaluatorIndependentFromFile.cc
		      evaluator/EvaluatorIndependentFunction.cc
		      evaluator/EvaluatorIndependentTensorFunction.cc
		      evaluator/EvaluatorPrimary.cc
                      evaluator/EvaluatorSecondary.cc
                      evaluator/EvaluatorSecondaryMonotype.cc
                      evaluator/EvaluatorSecondaryMonotypeFromFunction.cc
                      evaluator/Evaluator_OperatorApply.cc
#                      evaluator/Evaluator_PDE_Diffusion.cc
#                      evaluator/Evaluator_PDE_Accumulation.cc
                      evaluator/Evaluator_BCs.cc
                      evaluator/EvaluatorCellVolume.cc
#		      evaluator/cell_volume_evaluator.cc
#		      evaluator/deforming_cell_volume_evaluator.cc
#		      evaluator/rank_evaluator.cc
		      io/Debugger.cc
		      io/Checkpoint.cc
		      io/Visualization.cc
#		      io/Observable.cc
#		      io/UnstructuredObservations.cc
		      State.cc
                   HEADERS ${state_inc_files}
                   LINK_LIBS atk functions output mesh whetstone operators ${state_tpl_libs})
set_property(GLOBAL PROPERTY STATE_SOURCE_DIR ${MPC_SOURCE_DIR})
set_property(GLOBAL PROPERTY STATE_BINARY_DIR ${MPC_BINARY_DIR})

#
# Install Target
#
add_install_include_file(${state_inc_files})

register_evaluator_with_factory(HEADERFILE evaluator/EvaluatorPrimary_reg.hh LISTNAME REGISTER_AMANZI_STATE_EVALUATORS)
register_evaluator_with_factory(HEADERFILE evaluator/EvaluatorIndependentFunction_reg.hh LISTNAME REGISTER_AMANZI_STATE_EVALUATORS)
register_evaluator_with_factory(HEADERFILE evaluator/EvaluatorIndependentTensorFunction_reg.hh LISTNAME REGISTER_AMANZI_STATE_EVALUATORS)
register_evaluator_with_factory(HEADERFILE evaluator/EvaluatorIndependentFromFile_reg.hh LISTNAME REGISTER_AMANZI_STATE_EVALUATORS)
register_evaluator_with_factory(HEADERFILE evaluator/Evaluator_OperatorApply_reg.hh LISTNAME REGISTER_AMANZI_STATE_EVALUATORS)
#register_evaluator_with_factory(HEADERFILE evaluator/Evaluator_PDE_Diffusion_reg.hh LISTNAME REGISTER_AMANZI_STATE_EVALUATORS)
#register_evaluator_with_factory(HEADERFILE evaluator/Evaluator_PDE_Accumulation_reg.hh LISTNAME REGISTER_AMANZI_STATE_EVALUATORS)
register_evaluator_with_factory(HEADERFILE evaluator/Evaluator_BCs_reg.hh LISTNAME REGISTER_AMANZI_STATE_EVALUATORS)
register_evaluator_with_factory(HEADERFILE evaluator/Evaluator_BoundaryFunction_reg.hh LISTNAME REGISTER_AMANZI_STATE_EVALUATORS)
register_evaluator_with_factory(HEADERFILE evaluator/EvaluatorSecondaryMonotypeAdditive_reg.hh LISTNAME REGISTER_AMANZI_STATE_EVALUATORS)
register_evaluator_with_factory(HEADERFILE evaluator/EvaluatorSecondaryMonotypeMultiplicative_reg.hh LISTNAME REGISTER_AMANZI_STATE_EVALUATORS)
register_evaluator_with_factory(HEADERFILE evaluator/EvaluatorCellVolume_reg.hh LISTNAME REGISTER_AMANZI_STATE_EVALUATORS)
#register_evaluator_with_factory(HEADERFILE evaluator/EvaluatorSecondaryMonotypeFromFunction_reg.hh LISTNAME REGISTER_AMANZI_STATE_EVALUATORS)
# register_evaluator_with_factory(HEADERFILE evaluator/rank_evaluator_reg.hh          LISTNAME REGISTER_AMANZI_STATE_EVALUATORS)
# register_evaluator_with_factory(HEADERFILE evaluator/deforming_cell_volume_evaluator_reg.hh LISTNAME REGISTER_AMANZI_STATE_EVALUATORS)

generate_evaluators_registration_header(HEADERFILE state_evaluators_registration.hh LISTNAME REGISTER_AMANZI_STATE_EVALUATORS  INSTALL True)


if (BUILD_TESTS)
    # Add UnitTest includes
    include_directories(${UnitTest_INCLUDE_DIRS})

    # Add State include directory. 
    include_directories(${STATE_SOURCE_DIR})

    # Mesh factory includes
    include_directories(${MESH_FACTORY_SOURCE_DIR})

    # Copy test subdirectory for out of source builds
    if (NOT ("${STATE_SOURCE_DIR}" STREQUAL "${STATE_BINARY_DIR}"))
        execute_process(COMMAND ${CMAKE_COMMAND} -E 
          copy_directory ${STATE_SOURCE_DIR}/test ${STATE_BINARY_DIR}/test) 
    endif()

    add_amanzi_test(data_factory data_factory
      KIND unit
      SOURCE test/Main.cc test/data_factory.cc
      LINK_LIBS error_handling atk mesh operators data_structures mesh_functions state ${state_tpl_list} ${UnitTest_LIBRARIES})

    add_amanzi_test(record record
      KIND unit
      SOURCE test/Main.cc test/record.cc
      LINK_LIBS atk mesh operators data_structures mesh_functions state ${state_tpl_list} ${UnitTest_LIBRARIES})

    add_amanzi_test(state_mesh state_mesh
      KIND unit
      SOURCE test/Main.cc test/state_mesh.cc
      LINK_LIBS atk mesh_factory mesh operators data_structures mesh_functions state ${state_tpl_list} ${UnitTest_LIBRARIES})

    add_amanzi_test(state_as_container state_as_container
      KIND unit
      SOURCE test/Main.cc test/state_container.cc
      LINK_LIBS atk mesh_factory mesh operators data_structures mesh_functions operators state ${state_tpl_list} ${UnitTest_LIBRARIES})

    add_amanzi_test(state_vis state_vis
      KIND int
      SOURCE test/Main.cc test/state_vis.cc
      LINK_LIBS state mesh_factory mesh operators data_structures mesh_functions
                functions ${state_tpl_list} ${UnitTest_LIBRARIES})

    add_amanzi_test(state_restart state_restart
      KIND int
      SOURCE test/Main.cc test/state_restart.cc
      LINK_LIBS state mesh_factory mesh operators data_structures mesh_functions
                functions ${state_tpl_list} ${UnitTest_LIBRARIES})
              
    add_amanzi_test(state_init state_init
      KIND int
      SOURCE test/Main.cc test/state_init.cc
      LINK_LIBS state mesh_factory mesh operators data_structures mesh_functions functions ${state_tpl_list} ${UnitTest_LIBRARIES})

    add_amanzi_test(state_extensibility state_extensibility
      KIND int
      SOURCE test/Main.cc test/state_extensibility.cc
      LINK_LIBS state mesh_factory mesh operators data_structures mesh_functions functions ${state_tpl_list} ${UnitTest_LIBRARIES})
    

    add_amanzi_test(state_evaluators state_evaluators
      KIND int
      SOURCE test/Main.cc test/state_evaluators.cc
      LINK_LIBS state mesh_factory mesh operators data_structures mesh_functions functions ${UnitTest_LIBRARIES})

    add_amanzi_test(state_evaluators_cv state_evaluators_cv
      KIND int
      SOURCE test/Main.cc test/state_evaluators_cv.cc 
      LINK_LIBS state mesh_factory mesh operators data_structures mesh_functions functions ${UnitTest_LIBRARIES})

    add_amanzi_test(state_dag state_dag
      KIND int
      SOURCE test/Main.cc test/state_dag.cc
      LINK_LIBS state mesh_factory operators mesh data_structures mesh_functions functions ${UnitTest_LIBRARIES})

    add_amanzi_test(state_evaluators_op state_evaluators_op
      KIND int
      SOURCE test/Main.cc test/state_evaluators_op.cc
      LINK_LIBS state mesh_factory mesh operators data_structures mesh_functions functions operators ${UnitTest_LIBRARIES})

    add_amanzi_test(state_evaluators_op_inverse state_evaluators_op_inverse
      KIND int
      SOURCE test/Main.cc test/state_evaluators_op_inverse.cc
      LINK_LIBS state mesh_factory mesh data_structures mesh_functions functions ${UnitTest_LIBRARIES})


endif()

