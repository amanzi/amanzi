# -*- mode: cmake -*-

#
#  ATS
#    Flow PK class
#

# ATS include directories


# include_directories(${ATS_SOURCE_DIR}/src/data_structures)
# include_directories(${ATS_SOURCE_DIR}/src/state)
# include_directories(${ATS_SOURCE_DIR}/src/constitutive_relations/eos)
# include_directories(${ATS_SOURCE_DIR}/src/constitutive_relations/ewc)
# include_directories(${ATS_SOURCE_DIR}/src/constitutive_relations/generic_evaluators)
# include_directories(${ATS_SOURCE_DIR}/src/pks)
# include_directories(${ATS_SOURCE_DIR}/src/pks/mpc)
# include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations)
# include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/porosity)
# include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/elevation)
# include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/wrm)
# include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/overland_conductivity)
# include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/water_content)
# include_directories(${ATS_SOURCE_DIR}/src/pks/flow/richards)
# include_directories(${ATS_SOURCE_DIR}/src/pks/flow/permafrost)
# include_directories(${ATS_SOURCE_DIR}/src/pks/flow/overland)
# include_directories(${ATS_SOURCE_DIR}/src/pks/flow/overland_pressure)
#include_directories(${ATS_SOURCE_DIR}/src/pks/transport/transport_amanzi)


add_library(pk_morphology mpc_morphology_pk.cc)

install(TARGETS pk_morphology  DESTINATION lib)





