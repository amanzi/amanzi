# -*- mode: cmake -*-

#
#  ATS
#    Deformation PK class
#

#add_subdirectory(constitutive_relations)
#add_subdirectory(prescribed_deformation)
#add_subdirectory(prescribed_volumetric_deformation)
add_subdirectory(volumetric_deformation)

#============================================================================
# register evaluators/facotries/pks

# constitutive_relations/porosity/

register_evaluator_with_factory(
  HEADERFILE constitutive_relations/porosity/porosity_evaluator_reg.hh
  LISTNAME  ATS_POROSITY_EVAL)

generate_evaluators_registration_header(
  HEADERFILE deform_constitutive_relations_porosity_registration.hh
  LISTNAME   ATS_POROSITY_EVAL
  INSTALL    True
  )


# prescribed_deformation/

# register_evaluator_with_factory(
#   HEADERFILE prescribed_deformation/prescribed_deformation_reg.hh
#   LISTNAME  PRESCRIBED_DEFORMATION_EVAL)

# generate_evaluators_registration_header(
#   HEADERFILE deform_prescribed_deformation_registration.hh
#   LISTNAME   PRESCRIBED_DEFORMATION_EVAL
#   INSTALL    True
#   )

# volumetric_deformation/

register_evaluator_with_factory(
  HEADERFILE volumetric_deformation/volumetric_deformation_reg.hh
  LISTNAME  VOLUMETRIC_DEFORMATION_EVAL)

generate_evaluators_registration_header(
  HEADERFILE deform_volumetric_deformation_registration.hh
  LISTNAME   VOLUMETRIC_DEFORMATION_EVAL
  INSTALL    True
  )

  