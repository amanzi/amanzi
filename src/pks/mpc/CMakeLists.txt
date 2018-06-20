# -*- mode: cmake -*-

#
#  ATS
#    Flow PK class
#

# ATS include directories


include_directories(${ATS_SOURCE_DIR}/src/data_structures)
include_directories(${ATS_SOURCE_DIR}/src/state)
include_directories(${ATS_SOURCE_DIR}/src/constitutive_relations/eos)
include_directories(${ATS_SOURCE_DIR}/src/constitutive_relations/ewc)
include_directories(${ATS_SOURCE_DIR}/src/constitutive_relations/generic_evaluators)
include_directories(${ATS_SOURCE_DIR}/src/pks)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/porosity)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/elevation)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/wrm)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/overland_conductivity)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/water_content)

include_directories(${ATS_SOURCE_DIR}/src/pks/flow/richards)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/permafrost)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/overland)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/overland_pressure)
include_directories(${ATS_SOURCE_DIR}/src/pks/energy/constitutive_relations/internal_energy)
include_directories(${ATS_SOURCE_DIR}/src/pks/energy/base)
include_directories(${ATS_SOURCE_DIR}/src/pks/energy/two_phase)
include_directories(${ATS_SOURCE_DIR}/src/pks/transport/transport_amanzi)
include_directories(${ATS_SOURCE_DIR}/src/operators/advection)
include_directories(${ATS_SOURCE_DIR}/src/operators/divgrad/upwind_scheme)

add_library(mpc
  weak_mpc.cc
  operator_split_mpc.cc
  weak_mpc_semi_coupled.cc
  weak_mpc_semi_coupled_deform.cc
  mpc_coupled_cells.cc
  pk_mpcsubcycled_ats.cc
#  mpc_flowreactivetransport_pk.cc
#  mpc_reactivetransport_pk.cc
  mpc_surface_subsurface_helpers.cc
  weak_mpc_semi_coupled_helper.cc
  mpc_weak_subgrid.cc
  mpc_delegate_ewc.cc
##  mpc_delegate_ewc_surface.cc
  mpc_delegate_ewc_subsurface.cc
  mpc_delegate_water.cc
  mpc_coupled_water.cc
  mpc_coupled_water_split_flux.cc
#  mpc_coupled_water_neumann_splitting.cc
#  mpc_coupled_water_dirichlet_splitting.cc
  mpc_subsurface.cc
##  mpc_permafrost3.cc
  mpc_permafrost4.cc
##  mpc_coupled_transport.cc
)

install(TARGETS mpc DESTINATION lib)


# register factories
register_evaluator_with_factory(
  HEADERFILE weak_mpc_reg.hh
  LISTNAME   MPC_PK_FAC
  )

register_evaluator_with_factory(
  HEADERFILE strong_mpc_reg.hh
  LISTNAME   MPC_PK_FAC
  )

register_evaluator_with_factory(
  HEADERFILE mpc_predictor_corrector_reg.hh
  LISTNAME   MPC_PK_FAC
  )

register_evaluator_with_factory(
  HEADERFILE operator_split_mpc_reg.hh
  LISTNAME   MPC_PK_FAC
  )

register_evaluator_with_factory(
  HEADERFILE pk_mpcsubcycled_ats_reg.hh
  LISTNAME   MPC_PK_FAC
  )

# register_evaluator_with_factory(
#   HEADERFILE mpc_flowreactivetransport_pk_reg.hh
#   LISTNAME   MPC_PK_FAC
#   )

# register_evaluator_with_factory(
#   HEADERFILE mpc_reactivetransport_pk_reg.hh
#   LISTNAME   MPC_PK_FAC
#   )

    
register_evaluator_with_factory(
  HEADERFILE mpc_subsurface_reg.hh
  LISTNAME   MPC_PK_FAC
  )

register_evaluator_with_factory(
  HEADERFILE mpc_coupled_cells_reg.hh
  LISTNAME   MPC_PK_FAC
  )

register_evaluator_with_factory(
  HEADERFILE mpc_permafrost_reg.hh
  LISTNAME   MPC_PK_FAC
  )

register_evaluator_with_factory(
  HEADERFILE mpc_coupled_water_reg.hh
  LISTNAME   MPC_PK_FAC
  )

register_evaluator_with_factory(
  HEADERFILE mpc_coupled_water_split_flux_reg.hh
  LISTNAME   MPC_PK_FAC
  )

# register_evaluator_with_factory(
#   HEADERFILE mpc_coupled_water_neumann_splitting_reg.hh
#   LISTNAME   MPC_PK_FAC
#   )

# register_evaluator_with_factory(
#   HEADERFILE mpc_coupled_water_dirichlet_splitting_reg.hh
#   LISTNAME   MPC_PK_FAC
#   )

# register_evaluator_with_factory(
#    HEADERFILE coupled_transport/mpc_coupled_transport_reg.hh
#    LISTNAME   MPC_PK_FAC
#   )


register_evaluator_with_factory(
  HEADERFILE weak_mpc_semi_coupled_reg.hh
  LISTNAME   MPC_PK_FAC
  )

register_evaluator_with_factory(
  HEADERFILE mpc_weak_subgrid_reg.hh
  LISTNAME   MPC_PK_FAC
  )

register_evaluator_with_factory(
  HEADERFILE weak_mpc_semi_coupled_deform_reg.hh
  LISTNAME   MPC_PK_FAC
  )

generate_evaluators_registration_header(
  HEADERFILE mpc_registration.hh
  LISTNAME   MPC_PK_FAC
  INSTALL    True
)


include_evaluators_directories(LISTNAME ATS_POROSITY_EVAL_INCLUDES)
include_evaluators_directories(LISTNAME CONSTITUTIVE_RELATIONS_EOS_EVALUATORS_INCLUDES)
include_evaluators_directories(LISTNAME CONSTITUTIVE_RELATIONS_SURFACE_SUBSURFACE_FLUXES_EVALUATORS_INCLUDES)
include_evaluators_directories(LISTNAME CONSTITUTIVE_RELATIONS_GENERIC_EVALUATORS_INCLUDES)
include_evaluators_directories(LISTNAME FLOW_CONSTITUTIVE_RELATIONS_POROSITY_EVALUATORS_INCLUDES)
include_evaluators_directories(LISTNAME WRM_MODELS_EVALUATORS_FACTORIES_INCLUDES)
include_evaluators_directories(LISTNAME WMR_EVALUATORS_FACTORIES_INCLUDES)
include_evaluators_directories(LISTNAME RICHARDS_EVAL_FAC_INCLUDES)
include_evaluators_directories(LISTNAME PERMAFROST_EVAL_FAC_INCLUDES)
include_evaluators_directories(LISTNAME ENERGY_CONST_REL_INTERNAL_ENERGY_REG_INCLUDES)
include_evaluators_directories(LISTNAME ENERGY_CONST_REL_SOURCE_TERMS_REG_INCLUDES)
include_evaluators_directories(LISTNAME ENERGY_CONST_REL_THERMAL_CONDUCTIVITY_REG_INCLUDES)
include_evaluators_directories(LISTNAME ENERGY_TWO_PHASE_REG_INCLUDES)
include_evaluators_directories(LISTNAME ENERGY_THREE_PHASE_REG_INCLUDES)
include_evaluators_directories(LISTNAME MPC_PK_FAC_INCLUDES)


#add_subdirectory(coupled_transport)

