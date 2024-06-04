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
include_directories(${OUTPUT_SOURCE_DIR})
include_directories(${SOLVERS_SOURCE_DIR})
include_directories(${STATE_SOURCE_DIR}/data)
include_directories(${STATE_SOURCE_DIR}/evaluators)
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
  StateDefs.hh
  Debugger.hh
  data/Data.hh
  data/Data_Helpers.hh
  data/Data_Impl.hh
  data/DataFactory.hh
  data/DataFactory_Impl.hh
  data/Operator_DataHelpers.hh
  data/Record.hh
  data/RecordSet.hh
  evaluators/Evaluator.hh
  evaluators/Evaluator_Factory.hh
  evaluators/EvaluatorIndependent.hh
  evaluators/EvaluatorIndependentFromFile.hh
  evaluators/EvaluatorIndependentFunction.hh
  evaluators/EvaluatorIndependentConstant.hh
  evaluators/EvaluatorIndependentTensorFunction.hh
  evaluators/EvaluatorIndependentPatchFunction.hh
  evaluators/EvaluatorPrimary.hh
  evaluators/EvaluatorPrimaryStaticMesh.hh
  evaluators/EvaluatorSecondary.hh
  evaluators/EvaluatorSecondaryMonotype.hh
  evaluators/EvaluatorSecondaryMonotypeAdditive.hh
  evaluators/EvaluatorSecondaryMonotypeMultiplicative.hh
  evaluators/EvaluatorSecondaryMonotypeFromFunction.hh
  evaluators/EvaluatorSecondaryMeshedQuantity.hh
  evaluators/EvaluatorModelLauncher.hh
  evaluators/EvaluatorModelCV.hh
  evaluators/EvaluatorModelCVByMaterial.hh
  evaluators/EvaluatorMultiDOFModelCVByMaterial.hh
  evaluators/EvaluatorModelPatch.hh
  evaluators/EvaluatorAggregateBCs.hh
  evaluators/EvaluatorSecondaryVectorAsPatch.hh
  State.hh
  Observable.hh
  UnstructuredObservations.hh
  Visualization.hh
  VisualizationDomainSet.hh
  Checkpoint.hh
  # chemistry_data.hh 
  state_evaluators_registration.hh
  )

set(state_src_files
  Debugger.cc
  data/Record.cc
  data/RecordSet.cc
  data/Data_Helpers.cc
  evaluators/Evaluator.cc
  evaluators/Evaluator_Factory.cc
  evaluators/EvaluatorIndependent.cc
  evaluators/EvaluatorIndependentFromFile.cc
  evaluators/EvaluatorIndependentFunction.cc
  evaluators/EvaluatorIndependentConstant.cc
  evaluators/EvaluatorIndependentTensorFunction.cc
  evaluators/EvaluatorIndependentPatchFunction.cc
  evaluators/EvaluatorPrimary.cc
  evaluators/EvaluatorSecondary.cc
  evaluators/EvaluatorSecondaryMonotype.cc
  evaluators/EvaluatorSecondaryMonotypeFromFunction.cc
  evaluators/EvaluatorAggregateBCs.cc
  evaluators/EvaluatorSecondaryVectorAsPatch.cc
  State.cc
  Observable.cc
  UnstructuredObservations.cc
  Visualization.cc
  VisualizationDomainSet.cc
  Checkpoint.cc
  # chemistry_data.cc
  )  

set(state_link_libs  data_structures functions mesh_functions output mesh atk ${Boost_LIBRARIES})

add_amanzi_library(state 
  SOURCE ${state_src_files}
  HEADERS ${state_inc_files}
  LINK_LIBS ${state_link_libs})

# set_property(GLOBAL PROPERTY STATE_SOURCE_DIR ${MPC_SOURCE_DIR})
# set_property(GLOBAL PROPERTY STATE_BINARY_DIR ${MPC_BINARY_DIR})
if (APPLE AND BUILD_SHARED_LIBS)
  set_target_properties(state PROPERTIES LINK_FLAGS "-Wl,-undefined,dynamic_lookup")
endif()			  

#
# Install Target
#
add_install_include_file(${state_inc_files})

    # Add State include directory. 
    include_directories(${STATE_SOURCE_DIR})

    # Mesh factory includes
    include_directories(${MESH_EXTRACTED_SOURCE_DIR})
    include_directories(${MESH_LOGICAL_SOURCE_DIR})
    include_directories(${MESH_FACTORY_SOURCE_DIR})

if (BUILD_TESTS)

    # Add UnitTest includes
    include_directories(${UnitTest_INCLUDE_DIRS})
    set(state_test_link_libs state mesh_factory mesh data_structures mesh_functions functions
                              ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})


    # Copy test subdirectory for out of source builds
    if (NOT ("${STATE_SOURCE_DIR}" STREQUAL "${STATE_BINARY_DIR}"))
      file(GLOB DataFiles "${STATE_SOURCE_DIR}/test/*.xml"
                          "${STATE_SOURCE_DIR}/test/*.exo"
                          "${STATE_SOURCE_DIR}/test/*.gold")
      file(COPY ${DataFiles} DESTINATION ${STATE_BINARY_DIR}/test/)
    endif()

    # tests the "any" inside data objects
    add_amanzi_test(state_data_factory state_data_factory
                    KIND unit
                    SOURCE test/Main.cc 
                           test/state_data_factory.cc
                    LINK_LIBS ${state_test_link_libs})

    # tests record, recordset, and the underlying data objects
    add_amanzi_test(state_record state_record
                    KIND unit
                    SOURCE test/Main.cc
                           test/state_record.cc
                    LINK_LIBS ${state_test_link_libs})

    # tests state as a container of data                 
    add_amanzi_test(state_container state_container
                    KIND unit
                    SOURCE test/Main.cc
                           test/state_container.cc
                           test/state_extensibility.cc
                    LINK_LIBS ${state_test_link_libs})

    # tests state as a container of meshes
    add_amanzi_test(state_mesh state_mesh
                    KIND unit
                    SOURCE test/Main.cc
                           test/state_mesh.cc
                    LINK_LIBS ${state_test_link_libs})

    # tests the ability of state to initialize data                  
    add_amanzi_test(state_init state_init
                    KIND int
                    SOURCE test/Main.cc
         		   test/state_init.cc
                    LINK_LIBS ${state_test_link_libs})

    # tests the correctness of operator-assigning data by tag
    add_amanzi_test(state_copy state_copy
                    KIND int
                    SOURCE test/Main.cc test/state_copy.cc
                    LINK_LIBS ${state_test_link_libs})

    # tests observations, observation evaluators, and parallel observations
    add_amanzi_test(state_observations state_observations
                    KIND int
                    SOURCE test/Main.cc
         		   test/state_observations.cc
                    LINK_LIBS ${state_test_link_libs})

    if(ENABLE_MESH_MSTK)
      add_amanzi_test(state_observations_np2 state_observations NPROCS 2 KIND unit)
    endif()

    # tests the creation and writing of vis files
    add_amanzi_test(state_vis state_vis
                    KIND int
                    SOURCE test/Main.cc
                           test/state_vis.cc
                    LINK_LIBS ${state_test_link_libs})
                         
    # tests the creation, writing, and reading of checkpoint files on MPI_COMM_WORLD
    add_amanzi_test(state_restart state_restart
                    KIND int
                    SOURCE
          	          test/Main.cc
         		  test/state_restart.cc
                    LINK_LIBS ${state_test_link_libs})

    if(ENABLE_MESH_MSTK)
      # tests the creation, writing, and reading of checkpoint files on other comms
      add_amanzi_test(state_restart2 state_restart2
                    KIND int
                    SOURCE
          	          test/Main.cc
         		  test/state_restart2.cc
                    NPROCS 2
                    LINK_LIBS ${state_test_link_libs})
    endif() 

    # tests basic functionality of independent, primary, and secondary var evals
    add_amanzi_test(state_evaluators state_evaluators
                    KIND int
                    SOURCE test/Main.cc
                           test/state_evaluators.cc
                    LINK_LIBS ${state_test_link_libs})

    # same as the above, with CompositeVectors and additional vector functionality
    add_amanzi_test(state_evaluators_cv state_evaluators_cv
                    KIND int
                    SOURCE test/Main.cc
                           test/state_evaluators_cv.cc
                    LINK_LIBS ${state_test_link_libs})

    # tests Evaluators IndependentFromFunction and SecondaryFromFunction
    add_amanzi_test(state_functional_evaluators state_functional_evaluators
                    KIND int
                    SOURCE test/Main.cc
         	           test/state_functional_evaluators.cc
                    LINK_LIBS ${state_test_link_libs})

    # tests basic construction of and use of dag
    add_amanzi_test(state_dag state_dag
                    KIND int
                    SOURCE test/Main.cc test/state_dag.cc
                    LINK_LIBS ${state_test_link_libs})

    # tests Evaluators SecondaryMonotype<CV> in a dag
    add_amanzi_test(state_dag_cv state_dag_cv
                    KIND int
                    SOURCE test/Main.cc
         	           test/state_dag_cv.cc
                    LINK_LIBS ${state_test_link_libs})

    # tests Evaluators Model, ModelByMaterial in a dag
    add_amanzi_test(state_dag_models state_dag_models
                    KIND int
                    SOURCE test/Main.cc
         	           test/state_dag_models.cc
                    LINK_LIBS ${state_test_link_libs})

    # tests Evaluators on patches
    add_amanzi_test(state_domain_function_evaluators state_domain_function_evaluators
                    KIND int
                    SOURCE test/Main.cc
         	           test/state_domain_function_evaluators.cc
                    LINK_LIBS ${state_test_link_libs})

endif()

