# -*- mode: cmake -*-

#
#  Amanzi
#   Advection operators
#
# advection operators, including donor upwind, 2nd order with limiters, etc

add_library(advection advection.cc advection_donor_upwind.cc advection_factory.cc)
