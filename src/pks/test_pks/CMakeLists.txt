# -*- mode: cmake -*-

#
#  ATS
#    Flow PK class
#

add_subdirectory(divgrad_test)


#================================================
# register evaluators/factories/pks

# divgrad_test

register_evaluator_with_factory(
  HEADERFILE divgrad_test/divgrad_test_reg.hh
  LISTNAME   TEST_PKS_DIVGRAD_TEST_REG
  )

generate_evaluators_registration_header(
  HEADERFILE test_pks_divgrad_test_registration.hh
  LISTNAME   TEST_PKS_DIVGRAD_TEST_REG
  INSTALL    True
  )

