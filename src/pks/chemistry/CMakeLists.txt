# -*- mode: cmake -*-

#
#  ATS
#    Chemistry Amanzi
#

# ATS include directories

add_subdirectory(chemistry_amanzi)

# chemistry_amanzi

register_evaluator_with_factory(
   HEADERFILE chemistry_amanzi/Amanzi_PK_reg.hh
   LISTNAME   CHEMISTRY_AMANZI_REG
   )

generate_evaluators_registration_header(
   HEADERFILE chemistry_amanzi_registration.hh
   LISTNAME   CHEMISTRY_AMANZI_REG
   INSTALL    True
   )