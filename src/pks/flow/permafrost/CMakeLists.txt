# -*- mode: cmake -*-
include_directories(${ATS_SOURCE_DIR}/src/pks)
include_directories(${ATS_SOURCE_DIR}/src/factory)
include_directories(${ATS_SOURCE_DIR}/src/pks/bc_factory)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/bc_factory)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/eos)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/wrm)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/wrm/models)
include_directories(${ATS_SOURCE_DIR}/src/operators/divgrad)
include_directories(${ATS_SOURCE_DIR}/src/operators/divgrad/upwind_scheme)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/richards)

add_library(pk_flow_permafrost permafrost_pk.cc ) #interfrost.cc )

install(TARGETS pk_flow_permafrost DESTINATION lib)

# if ( BUILD_TESTS )
#     # Add UnitTest includes
#     set(amanzi_libs geometry mesh mesh_audit simple_mesh state whetstone mesh_functions)

#     include_directories(${Amanzi_TPL_UnitTest_INCLUDE_DIRS})

#     include_directories(${ATS_SOURCE_DIR}/src/pks/flow/permafrost)

#     # Test: miscaleneous transport routines
#     add_executable(permafrost_test test/Main.cc test/flow_test_class.cc test/permafrost_advance_simple.cc)
#     target_link_libraries(permafrost_test ${pk_factory} pk_flow_permafrost flow_bc_factory ${flow_relations_eos} ${flow_relations_wrm} divgrad amanzi_solvers amanzi_time_integration amanzi_whetstone amanzi_mesh amanzi_geometry amanzi_mesh_functions amanzi_mesh_factory amanzi_input amanzi_output amanzi_data_structures amanzi_state amanzi_data_structures ${Amanzi_TPL_UnitTest_LIBRARIES})

# endif()

