# -*- mode: cmake -*-

# ATS Surface balance PKs describe Evaporation, energy fluxes from
#  long/showtwave radiation, precip, etc etc etc


include_directories(${ATS_SOURCE_DIR}/src/constitutive_relations/surface_subsurface_fluxes)
include_directories(${ATS_SOURCE_DIR}/src/pks)
include_directories(constitutive_relations/SEB)

add_library(pk_surface_balance_SEB
       constitutive_relations/SEB/SnowEnergyBalance.cc
       constitutive_relations/SEB/SnowEnergyBalance_VPL.cc
       constitutive_relations/SEB/seb_physics_defs.cc
       constitutive_relations/SEB/seb_physics_funcs.cc
       SEB/surface_balance_SEB.cc
       SEB/surface_balance_evaluator.cc
       SEB_VPL/surface_balance_SEB_VPL.cc
       SEB_VPL/surface_balance_evaluator_VPL.cc
       implicit/surface_balance_implicit.cc)


install(TARGETS pk_surface_balance_SEB DESTINATION lib)

#================================================
# register evaluators/factories/pks

# SEB implicit
register_evaluator_with_factory(
  HEADERFILE implicit/surface_balance_implicit_reg.hh
  LISTNAME SURFACE_BALANCE_SEB_REG
)

# SEB
register_evaluator_with_factory(
  HEADERFILE SEB/surface_balance_SEB_reg.hh
  LISTNAME   SURFACE_BALANCE_SEB_REG
  )

# SEB VPL
register_evaluator_with_factory(
  HEADERFILE SEB_VPL/surface_balance_SEB_VPL_reg.hh
  LISTNAME   SURFACE_BALANCE_SEB_REG
  )

generate_evaluators_registration_header(
  HEADERFILE surface_balance_SEB_registration.hh
  LISTNAME   SURFACE_BALANCE_SEB_REG
  INSTALL    True
  )


if (BUILD_TESTS)
  include_directories(${Amanzi_TPL_Trilinos_INCLUDE_DIRS})
  include_directories(${Amanzi_TPL_UnitTest_INCLUDE_DIRS})
  include_directories(${Amanzi_TPL_Boost_INCLUDE_DIRS})
  include_directories(${ATS_SOURCE_DIR}/src/pks/surface_balance/constitutive_relations/SEB)

  add_executable(test_SEBNew_continuity
    constitutive_relations/SEB/test/test_SEBNew.cc
    constitutive_relations/SEB/test/main.cc)
  target_link_libraries(test_SEBNew_continuity
    pk_surface_balance_SEB
    amanzi_error_handling
    ${Amanzi_TPL_Teuchos_LIBRARIES}
    ${Amanzi_TPL_UnitTest_LIBRARIES}
    ${Amanzi_TPL_Boost_LIBRARIES})

  add_executable(test_SEBVPL_continuity
    constitutive_relations/SEB/test/test_SEB_VPL.cc
    constitutive_relations/SEB/test/main.cc)
  target_link_libraries(test_SEBVPL_continuity
    pk_surface_balance_SEB
    amanzi_error_handling
    ${Amanzi_TPL_Teuchos_LIBRARIES}
    ${Amanzi_TPL_UnitTest_LIBRARIES}
    ${Amanzi_TPL_Boost_LIBRARIES})

endif()