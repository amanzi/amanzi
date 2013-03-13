# -*- mode: cmake -*-

#
#  ATS
#    WRM for flow
#

include_directories(${ATS_SOURCE_DIR}/src/factory)

add_library(flow_relations_wrm_models
  wrm_factory.cc
  wrm_linear_system.cc
  wrm_van_genuchten.cc
  wrm_zero_relperm.cc
  pc_ice_water.cc
  wrm_permafrost_factory.cc
  wrm_old_permafrost_model.cc
  wrm_implicit_permafrost_model.cc
  wrm_implicit_permafrost_model2.cc
  )

if ( BUILD_TESTS )
    include_directories(${Amanzi_TPL_UnitTest_INCLUDE_DIRS})
    include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/wrm/models)

    add_executable(test_wrm test/main.cc test/test_vanGenuchten.cc test/test_old_permafrost.cc test/test_implicit_permafrost.cc test/test_implicit_permafrost2.cc)
    target_link_libraries(test_wrm flow_relations_wrm_models amanzi_error_handling ${Amanzi_TPL_UnitTest_LIBRARIES} ${Amanzi_TPL_Trilinos_LIBRARIES})
endif()



