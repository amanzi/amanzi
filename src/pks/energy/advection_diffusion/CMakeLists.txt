# -*- mode: cmake -*-

include_directories(${ATS_SOURCE_DIR}/src/pks)
include_directories(${ATS_SOURCE_DIR}/src/pks/energy/bc_factory)
include_directories(${ATS_SOURCE_DIR}/src/operators/advection)
include_directories(${ATS_SOURCE_DIR}/src/operators/divgrad)
add_library(pk_energy_advection_diffusion advection_diffusion_pk.cc advection_diffusion_ti.cc advection_diffusion_physics.cc)

if ( BUILD_TESTS )
  set(pk_factory -Wl,--whole-archive pk_factory -Wl,--no-whole-archive)

  # Add UnitTest includes
  set(amanzi_libs geometry mesh mesh_audit simple_mesh state whetstone mesh_functions)

  include_directories(${Amanzi_TPL_UnitTest_INCLUDE_DIRS})

  include_directories(${ATS_SOURCE_DIR}/src/pks/energy/advection_diffusion)

  # Test: miscaleneous transport routines
  add_executable(advection_diffusion_test test/Main.cc test/energy_test_class.cc test/advection_diffusion_advance_simple.cc)
  target_link_libraries(advection_diffusion_test ${pk_factory} amanzi_data_structures pk_energy_advection_diffusion energy_bc_factory advection divgrad amanzi_solvers amanzi_time_integration amanzi_whetstone amanzi_mesh amanzi_cell_geometry amanzi_mesh_functions amanzi_mesh_factory amanzi_input amanzi_output amanzi_data_structures amanzi_state ${Amanzi_TPL_UnitTest_LIBRARIES})

endif()

