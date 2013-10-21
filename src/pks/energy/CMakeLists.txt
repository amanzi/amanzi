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
