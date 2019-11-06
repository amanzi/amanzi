# -*- mode: cmake -*-

#
#  Generic Evaluators 
#   currently only header files
#

include_directories(${ATS_SOURCE_DIR}/src/factory)

add_library(generic_evaluators
    MultiplicativeEvaluator.cc
    AdditiveEvaluator.cc
    SubgridDisaggregateEvaluator.cc
    SubgridAggregateEvaluator.cc
    )

install(TARGETS generic_evaluators DESTINATION lib)
