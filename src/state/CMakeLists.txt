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
include_directories(${DBC_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${MESH_DATA_SOURCE_DIR})
include_directories(${MESH_FACTORY_SOURCE_DIR})
include_directories(${DATA_STRUCTURES_SOURCE_DIR})
include_directories(${OUTPUT_SOURCE_DIR})
include_directories(${FUNCS_SOURCE_DIR})
include_directories(${MFUNCS_SOURCE_DIR})
include_directories(${DBG_SOURCE_DIR})
include_directories(${DATA_STRUCTURES_SOURCE_DIR})
include_directories(${PLATFORM_INTERFACE_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${ASCEMIO_INCLUDE_DIR})


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
  cell_volume_evaluator.hh
  cell_volume_evaluator_reg.hh
  primary_variable_field_evaluator.hh
  primary_variable_field_evaluator_reg.hh
  secondary_variable_field_evaluator.hh
  secondary_variable_field_evaluator_reg.hh
  independent_variable_field_evaluator.hh
  independent_variable_field_evaluator_reg.hh
  deforming_cell_volume_evaluator.hh
  deforming_cell_volume_evaluator_reg.hh
  checkpoint.hh
  chemistry_data.hh
  Field_CompositeVector.hh
  Field_ConstantVector.hh
  FieldEvaluator_Factory.hh
  FieldEvaluator.hh
  Field.hh
  Field_Scalar.hh
  io_event.hh
  Observable.hh
  rank_evaluator.hh
  rank_evaluator_reg.hh
  StateDefs.hh
  State.hh
  UnstructuredObservations.hh
  visualization.hh)

set(state_tpl_list ${Boost_LIBRARIES} ${Teuchos_LIBRARIES} ${Epetra_LIBRARIES})
add_amanzi_library(state 
                   SOURCE State.cc
                          checkpoint.cc
                          visualization.cc
                          Observable.cc
                          UnstructuredObservations.cc
                          chemistry_data.cc
                          cell_volume_evaluator.cc
                          Field.cc
	                  Field_CompositeVector.cc
	                  Field_ConstantVector.cc
			  Field_Scalar.cc
			  FieldEvaluator.cc
	                  FieldEvaluator_Factory.cc
	                  io_event.cc Debugger.cc
                          rank_evaluator.cc
                          primary_variable_field_evaluator.cc
                          secondary_variable_field_evaluator.cc
                          secondary_variables_field_evaluator.cc
                          independent_variable_field_evaluator.cc
                          deforming_cell_volume_evaluator.cc
                   HEADERS ${state_inc_files}
                   LINK_LIBS atk functions output mesh ${state_tpl_libs})
set_property(GLOBAL PROPERTY STATE_SOURCE_DIR ${MPC_SOURCE_DIR})
set_property(GLOBAL PROPERTY STATE_BINARY_DIR ${MPC_BINARY_DIR})

#
# Install Target
#
add_install_include_file(${state_inc_files})

register_evaluator_with_factory(HEADERFILE cell_volume_evaluator_reg.hh   LISTNAME REGISTER_AMANZI_STATE_EVALUATORS)
register_evaluator_with_factory(HEADERFILE rank_evaluator_reg.hh          LISTNAME REGISTER_AMANZI_STATE_EVALUATORS)
register_evaluator_with_factory(HEADERFILE primary_variable_field_evaluator_reg.hh LISTNAME REGISTER_AMANZI_STATE_EVALUATORS)
register_evaluator_with_factory(HEADERFILE primary_variable_field_evaluator_reg.hh LISTNAME REGISTER_AMANZI_STATE_EVALUATORS)
register_evaluator_with_factory(HEADERFILE deforming_cell_volume_evaluator_reg.hh LISTNAME REGISTER_AMANZI_STATE_EVALUATORS)

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
    
    #set(amanzi_libs geometry mesh state)
    
    if ( ENABLE_STK_Mesh )
      add_amanzi_test(state_restart state_restart
                      KIND int
                      SOURCE
	  	          test/Main.cc
	 		  test/test_restart.cc
	 		  test/test_vis.cc
                      LINK_LIBS state mesh_factory mesh data_structures mesh_functions functions ${UnitTest_LIBRARIES})

      add_amanzi_test(state_copy state_copy
                      KIND int
                      SOURCE
	  	          test/Main.cc
	 		  test/test_copy.cc
                      LINK_LIBS state mesh_factory mesh data_structures mesh_functions functions ${UnitTest_LIBRARIES})
    endif()

endif()
