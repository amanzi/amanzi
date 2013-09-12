# -*- mode: cmake -*-
MACRO(WHOLE_ARCHIVE inlib)
  SET(${inlib} -Wl,--whole-archive ${inlib} -Wl,--no-whole-archive)
ENDMACRO(WHOLE_ARCHIVE inlib)

include_directories(${ATS_SOURCE_DIR}/src/pks)
include_directories(${ATS_SOURCE_DIR}/src/pks/mpc)
include_directories(${ATS_SOURCE_DIR}/src/factory)

include_directories(${Amanzi_TPL_MSTK_INCLUDE_DIRS})
add_definitions("-DMSTK_HAVE_MPI")


register_evaluator_with_factory(
  HEADERFILE porosity_evaluator_reg.hh 
  LISTNAME  ATS_POROSITY_EVAL)

add_library(deform_relations_porosity porosity_evaluator.cc)

