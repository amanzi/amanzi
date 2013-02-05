# -*- mode: cmake -*-
MACRO(WHOLE_ARCHIVE inlib)
  SET(${inlib} -Wl,--whole-archive ${inlib} -Wl,--no-whole-archive)
ENDMACRO(WHOLE_ARCHIVE inlib)

include_directories(${ATS_SOURCE_DIR}/src/pks)
include_directories(${ATS_SOURCE_DIR}/src/factory)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/bc_factory)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/overland_conductivity)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/elevation)
include_directories(${ATS_SOURCE_DIR}/src/operators/divgrad)
include_directories(${ATS_SOURCE_DIR}/src/operators/divgrad/upwind_scheme)

add_library(pk_flow_overland
        overland.cc overland_ti.cc overland_physics.cc
        surface_coupler_via_source_evaluator.cc)


# shouldn't this get picked up from ../constitutive_relations/eos/CMakeLists.txt?
WHOLE_ARCHIVE(flow_relations_eos)
WHOLE_ARCHIVE(flow_relations_wrm)
WHOLE_ARCHIVE(pk_factory)

# if ( BUILD_TESTS )
#     # Add UnitTest includes
#     include_directories(${Amanzi_TPL_MSTK_INCLUDE_DIRS})
#     add_definitions("-DMSTK_HAVE_MPI")

#     include_directories(${Amanzi_TPL_UnitTest_INCLUDE_DIRS})
#     include_directories(${ATS_SOURCE_DIR}/src/pks/flow/overland)

#     # Test: miscaleneous transport routines
#     add_executable(overland_test test/Main.cc test/flow_test_class.cc test/overland_advance_simple.cc)
#     target_link_libraries(overland_test ${pk_factory} pk_flow_overland flow_bc_factory amanzi_mesh_functions ${flow_relations_eos} ${flow_relations_wrm} divgrad amanzi_solvers amanzi_time_integration amanzi_whetstone amanzi_mesh amanzi_geometry amanzi_mesh_factory amanzi_input amanzi_output amanzi_data_structures amanzi_state amanzi_data_structures ${Amanzi_TPL_UnitTest_LIBRARIES})

# endif()

