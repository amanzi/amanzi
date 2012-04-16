# -*- mode: cmake -*-

include_directories(${ATS_SOURCE_DIR}/src/pks)
include_directories(${ATS_SOURCE_DIR}/src/operators/advection)
add_library(pk_transport_passive_tracer passive_tracer.cc passive_tracer_ti.cc passive_tracer_physics.cc)

if ( BUILD_TESTS )
  set(pk_factory -Wl,--whole-archive pk_factory -Wl,--no-whole-archive)

  # Add UnitTest includes
  set(amanzi_libs geometry mesh mesh_audit simple_mesh state transport whetstone mesh_functions)

  include_directories(${Amanzi_TPL_UnitTest_INCLUDE_DIRS})

  include_directories(${ATS_SOURCE_DIR}/src/pks/transport/passive_tracer)

  # Test: miscaleneous transport routines
  add_executable(transport_test test/Main.cc test/transport_test_class.cc test/transport_advance_simple.cc)
  target_link_libraries(transport_test ${pk_factory} amanzi_data_structures pk_transport_passive_tracer advection amanzi_solvers amanzi_time_integration amanzi_mesh amanzi_cell_geometry amanzi_mesh_functions amanzi_mesh_factory amanzi_input amanzi_output amanzi_data_structures amanzi_state ${Amanzi_TPL_UnitTest_LIBRARIES})

endif()
