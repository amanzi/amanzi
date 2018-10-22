# -*- mode: cmake -*-

# ATS Surface balance PKs describe Evaporation, energy fluxes from
#  long/showtwave radiation, precip, etc etc etc


include_directories(${ATS_SOURCE_DIR}/src/constitutive_relations/surface_subsurface_fluxes)
include_directories(${ATS_SOURCE_DIR}/src/pks)
include_directories(constitutive_relations/SEB)
include_directories(constitutive_relations/litter)
include_directories(base)

add_library(pk_surface_balance_SEB
       constitutive_relations/SEB/seb_physics_defs.cc
       constitutive_relations/SEB/seb_physics_funcs.cc
       constitutive_relations/SEB/albedo_evaluator.cc
       constitutive_relations/SEB/albedo_subgrid_evaluator.cc
       constitutive_relations/SEB/area_fractions_evaluator.cc
       constitutive_relations/SEB/area_fractions_subgrid_evaluator.cc
       constitutive_relations/SEB/seb_evaluator.cc
       constitutive_relations/SEB/seb_subgrid_evaluator.cc
       constitutive_relations/SEB/longwave_evaluator.cc
       constitutive_relations/litter/drainage_evaluator.cc
       constitutive_relations/litter/interception_evaluator.cc
       constitutive_relations/litter/latent_heat_evaluator.cc
       constitutive_relations/litter/latent_heat_model.cc       
       constitutive_relations/litter/evaporative_flux_relaxation_evaluator.cc       
       constitutive_relations/litter/evaporative_flux_relaxation_model.cc
       constitutive_relations/litter/micropore_macropore_flux_evaluator.cc       
       constitutive_relations/litter/micropore_macropore_flux_model.cc       
       constitutive_relations/litter/macropore_surface_flux_evaluator.cc       
       constitutive_relations/litter/macropore_surface_flux_model.cc       
       base/surface_balance_base.cc
       implicit_subgrid/surface_balance_implicit_subgrid.cc
)


install(TARGETS pk_surface_balance_SEB DESTINATION lib)

#================================================
# register evaluators/factories/pks

# SEB evalutors
register_evaluator_with_factory(
  HEADERFILE constitutive_relations/SEB/seb_evaluator_reg.hh
  LISTNAME SURFACE_BALANCE_SEB_REG
)

register_evaluator_with_factory(
  HEADERFILE constitutive_relations/SEB/seb_subgrid_evaluator_reg.hh
  LISTNAME SURFACE_BALANCE_SEB_REG
)

register_evaluator_with_factory(
  HEADERFILE constitutive_relations/SEB/albedo_evaluator_reg.hh
  LISTNAME SURFACE_BALANCE_SEB_REG
)
register_evaluator_with_factory(
  HEADERFILE constitutive_relations/SEB/albedo_subgrid_evaluator_reg.hh
  LISTNAME SURFACE_BALANCE_SEB_REG
)

register_evaluator_with_factory(
  HEADERFILE constitutive_relations/SEB/area_fractions_evaluator_reg.hh
  LISTNAME SURFACE_BALANCE_SEB_REG
)

register_evaluator_with_factory(
  HEADERFILE constitutive_relations/SEB/area_fractions_subgrid_evaluator_reg.hh
  LISTNAME SURFACE_BALANCE_SEB_REG
)

register_evaluator_with_factory(
  HEADERFILE constitutive_relations/SEB/longwave_evaluator_reg.hh
  LISTNAME SURFACE_BALANCE_SEB_REG
)

register_evaluator_with_factory(
  HEADERFILE constitutive_relations/litter/drainage_evaluator_reg.hh
  LISTNAME SURFACE_BALANCE_SEB_REG
)

register_evaluator_with_factory(
  HEADERFILE constitutive_relations/litter/interception_evaluator_reg.hh
  LISTNAME SURFACE_BALANCE_SEB_REG
)

register_evaluator_with_factory(
  HEADERFILE constitutive_relations/litter/latent_heat_evaluator_reg.hh
  LISTNAME SURFACE_BALANCE_SEB_REG
)

register_evaluator_with_factory(
  HEADERFILE constitutive_relations/litter/micropore_macropore_flux_evaluator_reg.hh
  LISTNAME SURFACE_BALANCE_SEB_REG
)

register_evaluator_with_factory(
  HEADERFILE constitutive_relations/litter/macropore_surface_flux_evaluator_reg.hh
  LISTNAME SURFACE_BALANCE_SEB_REG
)


register_evaluator_with_factory(
  HEADERFILE constitutive_relations/litter/evaporative_flux_relaxation_evaluator_reg.hh
  LISTNAME SURFACE_BALANCE_SEB_REG
)


register_evaluator_with_factory(
  HEADERFILE base/surface_balance_base_reg.hh
  LISTNAME SURFACE_BALANCE_SEB_REG
)

register_evaluator_with_factory(
  HEADERFILE implicit_subgrid/surface_balance_implicit_subgrid_reg.hh
  LISTNAME SURFACE_BALANCE_SEB_REG
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

  # Not actually tests -- these just write stuff to file for debugging --etc
  
  # add_executable(test_SEBNew_continuity
  #   constitutive_relations/SEB/test/test_SEBNew.cc
  #   constitutive_relations/SEB/test/main.cc)
  # target_link_libraries(test_SEBNew_continuity
  #   pk_surface_balance_SEB
  #   amanzi_error_handling
  #   ${Amanzi_TPL_Teuchos_LIBRARIES}
  #   ${Amanzi_TPL_UnitTest_LIBRARIES}
  #   ${Amanzi_TPL_Boost_LIBRARIES})

  # add_executable(test_SEBVPL_continuity
  #   constitutive_relations/SEB/test/test_SEB_VPL.cc
  #   constitutive_relations/SEB/test/main.cc)
  # target_link_libraries(test_SEBVPL_continuity
  #   pk_surface_balance_SEB
  #   amanzi_error_handling
  #   ${Amanzi_TPL_Teuchos_LIBRARIES}
  #   ${Amanzi_TPL_UnitTest_LIBRARIES}
  #   ${Amanzi_TPL_Boost_LIBRARIES})

endif()
