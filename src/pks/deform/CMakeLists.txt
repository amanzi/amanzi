# -*- mode: cmake -*-

#
#  ATS
#    Deformation PK class
#

add_subdirectory(volumetric_deformation)

#============================================================================
# register evaluators/facotries/pks


# volumetric_deformation/
register_evaluator_with_factory(
  HEADERFILE volumetric_deformation/volumetric_deformation_reg.hh
  LISTNAME  DEFORMATION_REG)


generate_evaluators_registration_header(
  HEADERFILE deformation_registration.hh
  LISTNAME   DEFORMATION_REG
  INSTALL    True
  )

  
