# -*- mode: cmake -*-

#
#  ATS
#

include_directories(${ATS_SOURCE_DIR}/src/constitutive_relations/eos)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/wrm)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/wrm/models)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/porosity)
include_directories(${ATS_SOURCE_DIR}/src/pks/energy/constitutive_relations/internal_energy)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/overland_conductivity)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/overland_pressure)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/icy_overland)

add_library(relations_ewc
#  thermal_richards_model.cc
  ewc_model_base.cc
  liquid_ice_model.cc
  permafrost_model.cc
  surface_ice_model.cc
  )


install(TARGETS relations_ewc DESTINATION lib )



