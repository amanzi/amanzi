# -*- mode: cmake -*-
include_directories(${ATS_SOURCE_DIR}/src/pks)
include_directories(${ATS_SOURCE_DIR}/src/pks/mpc)
include_directories(${ATS_SOURCE_DIR}/src/factory)

include_directories(${Amanzi_TPL_MSTK_INCLUDE_DIRS})
add_definitions("-DMSTK_HAVE_MPI")


add_library(deform_relations_porosity porosity_evaluator.cc)

install(TARGETS deform_relations_porosity DESTINATION lib)


