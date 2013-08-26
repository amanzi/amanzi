# -*- mode: cmake -*-
MACRO(WHOLE_ARCHIVE inlib)
  SET(${inlib} -Wl,--whole-archive ${inlib} -Wl,--no-whole-archive)
ENDMACRO(WHOLE_ARCHIVE inlib)

include_directories(${ATS_SOURCE_DIR}/src/pks)
include_directories(${ATS_SOURCE_DIR}/src/pks/mpc)
include_directories(${ATS_SOURCE_DIR}/src/factory)
include_directories(${ATS_SOURCE_DIR}/src/pks/deform/constitutive_relations/porosity)

include_directories(${Amanzi_TPL_MSTK_INCLUDE_DIRS})
add_definitions("-DMSTK_HAVE_MPI")

add_library(pk_prescribed_volumetric_deformation prescribed_volumetric_deformation.cc)

if ( BUILD_TESTS )

endif()
