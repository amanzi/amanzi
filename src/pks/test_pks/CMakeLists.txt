# -*- mode: cmake -*-

#
#  ATS
#    Flow PK class
#

# -*- mode: cmake -*-
include_directories(${ATS_SOURCE_DIR}/src/pks)
include_directories(${ATS_SOURCE_DIR}/src/factory)
include_directories(${ATS_SOURCE_DIR}/src/pks/bc_factory)
include_directories(${ATS_SOURCE_DIR}/src/operators/divgrad)

add_library(pk_test_pks divgrad_test/divgrad_test.cc
                        test_snow_dist/test_snow_dist.cc)

install(TARGETS pk_test_pks DESTINATION lib)


#================================================
# register evaluators/factories/pks

# divgrad_test

register_evaluator_with_factory(
  HEADERFILE divgrad_test/divgrad_test_reg.hh
  LISTNAME   TEST_PKS_REG
  )

register_evaluator_with_factory(
  HEADERFILE test_snow_dist/test_snow_dist_reg.hh
  LISTNAME   TEST_PKS_REG
  )

generate_evaluators_registration_header(
  HEADERFILE test_pks_registration.hh
  LISTNAME   TEST_PKS_REG
  INSTALL    True
  )

