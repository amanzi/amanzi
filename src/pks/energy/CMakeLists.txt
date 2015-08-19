# -*- mode: cmake -*-

#
#  ATS
#    Flow PK class
#

# ATS include directories

# additional common components
add_subdirectory(constitutive_relations)

# PKs
add_subdirectory(constant_temperature)
add_subdirectory(advection_diffusion)
add_subdirectory(base)
add_subdirectory(two_phase)
add_subdirectory(three_phase)
add_subdirectory(surface_ice)

#================================================
# register evaluators/factories/pks

# advection_diffusion/

register_evaluator_with_factory(
  HEADERFILE advection_diffusion/advection_diffusion_pk_reg.hh
  LISTNAME   ENERGY_ADVECTION_DIFFUSION_REG
  )

generate_evaluators_registration_header(
  HEADERFILE energy_advection_diffusion_registration.hh
  LISTNAME   ENERGY_ADVECTION_DIFFUSION_REG
  INSTALL    True
  )
 

# constant_temperature/

register_evaluator_with_factory(
  HEADERFILE constant_temperature/constant_temperature_reg.hh
  LISTNAME   ENERGY_CONSTANT_TEMPERATURE_REG
  )

generate_evaluators_registration_header(
  HEADERFILE energy_constant_temperature_registration.hh
  LISTNAME   ENERGY_CONSTANT_TEMPERATURE_REG
  INSTALL    True
  )

# constitutive_relations/source_terms/
register_evaluator_with_factory(
  HEADERFILE constitutive_relations/enthalpy/enthalpy_evaluator_reg.hh
  LISTNAME   ENERGY_CONST_REL_SOURCE_TERMS_REG
  )

# constitutive_relations/source_terms/
register_evaluator_with_factory(
  HEADERFILE constitutive_relations/source_terms/advected_energy_source_evaluator_reg.hh
  LISTNAME   ENERGY_CONST_REL_SOURCE_TERMS_REG
  )

generate_evaluators_registration_header(
  HEADERFILE energy_constitutive_relations_source_terms_registration.hh
  LISTNAME   ENERGY_CONST_REL_SOURCE_TERMS_REG
  INSTALL    True
  )

# constitutive_relations/internal_energy/

register_evaluator_with_factory(
  HEADERFILE constitutive_relations/internal_energy/iem_evaluator_reg.hh
  LISTNAME   ENERGY_CONST_REL_INTERNAL_ENERGY_REG
  )

register_evaluator_with_factory(
  HEADERFILE constitutive_relations/internal_energy/iem_factory_reg.hh
  LISTNAME   ENERGY_CONST_REL_INTERNAL_ENERGY_REG
  )

register_evaluator_with_factory(
  HEADERFILE constitutive_relations/internal_energy/iem_linear_reg.hh
  LISTNAME   ENERGY_CONST_REL_INTERNAL_ENERGY_REG
  )

register_evaluator_with_factory(
  HEADERFILE constitutive_relations/internal_energy/iem_quadratic_reg.hh
  LISTNAME   ENERGY_CONST_REL_INTERNAL_ENERGY_REG
  )

register_evaluator_with_factory(
  HEADERFILE constitutive_relations/internal_energy/iem_water_vapor_evaluator_reg.hh
  LISTNAME   ENERGY_CONST_REL_INTERNAL_ENERGY_REG
  )

generate_evaluators_registration_header(
  HEADERFILE energy_constitutive_relations_internal_energy_registration.hh
  LISTNAME   ENERGY_CONST_REL_INTERNAL_ENERGY_REG
  INSTALL    True
  )


# constitutive_relations/thermal_conductivity/

register_evaluator_with_factory(
  HEADERFILE constitutive_relations/thermal_conductivity/thermal_conductivity_twophase_factory_reg.hh
  LISTNAME   ENERGY_CONST_REL_THERMAL_CONDUCTIVITY_REG
  )

register_evaluator_with_factory(
  HEADERFILE constitutive_relations/thermal_conductivity/thermal_conductivity_threephase_factory_reg.hh
  LISTNAME   ENERGY_CONST_REL_THERMAL_CONDUCTIVITY_REG
  )

register_evaluator_with_factory(
  HEADERFILE constitutive_relations/thermal_conductivity/thermal_conductivity_threephase_peterslidard_reg.hh
  LISTNAME   ENERGY_CONST_REL_THERMAL_CONDUCTIVITY_REG
  )

register_evaluator_with_factory(
  HEADERFILE constitutive_relations/thermal_conductivity/thermal_conductivity_threephase_wetdry_reg.hh
  LISTNAME   ENERGY_CONST_REL_THERMAL_CONDUCTIVITY_REG
  )

register_evaluator_with_factory(
  HEADERFILE constitutive_relations/thermal_conductivity/thermal_conductivity_threephase_volume_averaged_reg.hh
  LISTNAME   ENERGY_CONST_REL_THERMAL_CONDUCTIVITY_REG
  )

register_evaluator_with_factory(
  HEADERFILE constitutive_relations/thermal_conductivity/thermal_conductivity_threephase_sutra_hacked_reg.hh
  LISTNAME   ENERGY_CONST_REL_THERMAL_CONDUCTIVITY_REG
  )

register_evaluator_with_factory(
  HEADERFILE constitutive_relations/thermal_conductivity/thermal_conductivity_twophase_peterslidard_reg.hh
  LISTNAME   ENERGY_CONST_REL_THERMAL_CONDUCTIVITY_REG
  )

register_evaluator_with_factory(
  HEADERFILE constitutive_relations/thermal_conductivity/thermal_conductivity_twophase_wetdry_reg.hh
  LISTNAME   ENERGY_CONST_REL_THERMAL_CONDUCTIVITY_REG
  )

generate_evaluators_registration_header(
  HEADERFILE energy_constitutive_relations_thermal_conductivity_registration.hh
  LISTNAME   ENERGY_CONST_REL_THERMAL_CONDUCTIVITY_REG
  INSTALL    True
  )

# surface_ice/

register_evaluator_with_factory(
  HEADERFILE surface_ice/energy_surface_ice_reg.hh
  LISTNAME   ENERGY_SURFACE_ICE_REG
  )

register_evaluator_with_factory(
  HEADERFILE surface_ice/surface_ice_energy_evaluator_reg.hh
  LISTNAME   ENERGY_SURFACE_ICE_REG
  )

generate_evaluators_registration_header(
  HEADERFILE energy_surface_ice_registration.hh
  LISTNAME   ENERGY_SURFACE_ICE_REG
  INSTALL    True
  )


# three_phase/

register_evaluator_with_factory(
  HEADERFILE three_phase/three_phase_reg.hh
  LISTNAME   ENERGY_THREE_PHASE_REG
  )
register_evaluator_with_factory(
  HEADERFILE three_phase/interfrost_energy_reg.hh
  LISTNAME   ENERGY_THREE_PHASE_REG
  )
register_evaluator_with_factory(
  HEADERFILE three_phase/interfrost_energy_evaluator_reg.hh
  LISTNAME   ENERGY_THREE_PHASE_REG
  )
register_evaluator_with_factory(
  HEADERFILE three_phase/three_phase_energy_evaluator_reg.hh
  LISTNAME   ENERGY_THREE_PHASE_REG
  )

generate_evaluators_registration_header(
  HEADERFILE energy_three_phase_registration.hh
  LISTNAME   ENERGY_THREE_PHASE_REG
  INSTALL    True
  )


# two_phase/

register_evaluator_with_factory(
  HEADERFILE two_phase/two_phase_reg.hh
  LISTNAME   ENERGY_TWO_PHASE_REG
  )

generate_evaluators_registration_header(
  HEADERFILE energy_two_phase_registration.hh
  LISTNAME   ENERGY_TWO_PHASE_REG
  INSTALL    True
  )
