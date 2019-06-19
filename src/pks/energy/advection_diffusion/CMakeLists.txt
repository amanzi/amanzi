# -*- mode: cmake -*-

include_directories(${ATS_SOURCE_DIR}/src/pks)
include_directories(${ATS_SOURCE_DIR}/src/pks/bc_factory)
include_directories(${ATS_SOURCE_DIR}/src/pks/energy/bc_factory)
include_directories(${ATS_SOURCE_DIR}/src/operators/advection)

add_library(pk_energy_advection_diffusion advection_diffusion_pk.cc advection_diffusion_ti.cc advection_diffusion_physics.cc)

install(TARGETS pk_energy_advection_diffusion DESTINATION lib)

#if ( BUILD_TESTS )
if (0)
  # Add UnitTest includes
  include_directories(${Amanzi_TPL_UnitTest_INCLUDE_DIRS})

  include_directories(${ATS_SOURCE_DIR}/src/pks/energy/advection_diffusion)

  # Test: miscaleneous transport routines
  add_executable(advection_diffusion_test test/Main.cc test/energy_test_class.cc test/advection_diffusion_advance_simple.cc)
  target_link_libraries(advection_diffusion_test 
    pk_energy_advection_diffusion
    bc_factory
    pk_bases
    amanzi_pks
    advection 
    amanzi_solvers amanzi_time_integration amanzi_operators
    amanzi_state amanzi_whetstone amanzi_mesh amanzi_geometry
    amanzi_functions amanzi_mesh_functions amanzi_mesh_factory
    amanzi_output amanzi_data_structures amanzi_atk ${Amanzi_TPL_UnitTest_LIBRARIES})

endif()

