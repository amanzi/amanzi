# -*- mode: cmake -*-

include_directories(${ATS_SOURCE_DIR}/src/pks/energy/bc_factory)
include_directories(${ATS_SOURCE_DIR}/src/pks/energy/constitutive_relations)
include_directories(${ATS_SOURCE_DIR}/src/operators/advection)
include_directories(${ATS_SOURCE_DIR}/src/operators/divgrad)
add_library(energy_air_water_rock air_water_rock_pk.cc air_water_rock_ti.cc air_water_rock_physics.cc)

if ( BUILD_TESTS )
    # Add UnitTest includes
    set(amanzi_libs geometry mesh mesh_audit simple_mesh state whetstone mesh_functions)

    include_directories(${Amanzi_TPL_UnitTest_INCLUDE_DIRS})

    include_directories(${ATS_SOURCE_DIR}/src/pks/energy/air_water_rock)

    # Test: miscaleneous transport routines
    add_executable(air_water_rock_test test/Main.cc test/energy_test_class.cc test/air_water_rock_advance_simple.cc)
    target_link_libraries(air_water_rock_test amanzi_data_structures energy_air_water_rock energy_bc_factory energy_constitutive_relations advection divgrad amanzi_solvers amanzi_time_integration amanzi_whetstone amanzi_mesh amanzi_cell_geometry amanzi_mesh_functions amanzi_mesh_factory amanzi_input amanzi_output amanzi_data_structures amanzi_state ${Amanzi_TPL_UnitTest_LIBRARIES})

endif()

