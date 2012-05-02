# -*- mode: cmake -*-
MACRO(WHOLE_ARCHIVE inlib)
  SET(${inlib} -Wl,--whole-archive ${inlib} -Wl,--no-whole-archive)
ENDMACRO(WHOLE_ARCHIVE inlib)

include_directories(${ATS_SOURCE_DIR}/src/pks)
include_directories(${ATS_SOURCE_DIR}/src/factory)
include_directories(${ATS_SOURCE_DIR}/src/pks/energy/bc_factory)
include_directories(${ATS_SOURCE_DIR}/src/pks/energy/constitutive_relations/thermal_conductivity)
include_directories(${ATS_SOURCE_DIR}/src/pks/energy/constitutive_relations/internal_energy)
include_directories(${ATS_SOURCE_DIR}/src/operators/advection)
include_directories(${ATS_SOURCE_DIR}/src/operators/divgrad)
add_library(pk_energy_three_phase three_phase_pk.cc three_phase_ti.cc three_phase_physics.cc)

# WHOLE_ARCHIVE(energy_relations_thermal_conductivity)

# if ( BUILD_TESTS )
#   set(pk_factory -Wl,--whole-archive pk_factory -Wl,--no-whole-archive)

#   include_directories(${Amanzi_TPL_UnitTest_INCLUDE_DIRS})

#   include_directories(${ATS_SOURCE_DIR}/src/pks/energy/three_phase)

#   # Test: miscaleneous transport routines
#   add_executable(three_phase_test test/Main.cc test/energy_test_class.cc test/three_phase_advance_simple.cc)
#   target_link_libraries(three_phase_test ${pk_factory} amanzi_data_structures pk_energy_three_phase energy_bc_factory ${energy_relations_thermal_conductivity} energy_relations_internal_energy advection divgrad amanzi_solvers amanzi_time_integration amanzi_whetstone amanzi_mesh amanzi_cell_geometry amanzi_mesh_functions amanzi_mesh_factory amanzi_input amanzi_output amanzi_data_structures amanzi_state ${Amanzi_TPL_UnitTest_LIBRARIES})

# endif()

