# -*- mode: cmake -*-

#
#  ATS
#    Deformation PK class
#

add_subdirectory(constitutive_relations)
add_subdirectory(prescribed_deformation)
#add_subdirectory(prescribed_volumetric_deformation)
add_subdirectory(volumetric_deformation)

generate_evaluators_registration_header(
  HEADERFILE deform_evaluators_reg.hh
  LISTNAME   ATS_POROSITY_EVAL)

