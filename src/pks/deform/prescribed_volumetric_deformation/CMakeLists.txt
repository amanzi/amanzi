# -*- mode: cmake -*-
include_directories(${ATS_SOURCE_DIR}/src/pks)
include_directories(${ATS_SOURCE_DIR}/src/pks/mpc)
include_directories(${ATS_SOURCE_DIR}/src/factory)
include_directories(${ATS_SOURCE_DIR}/src/pks/deform/constitutive_relations/porosity)

include_directories(${Amanzi_TPL_MSTK_INCLUDE_DIRS})
add_definitions("-DMSTK_HAVE_MPI")

add_library(pk_prescribed_volumetric_deformation prescribed_volumetric_deformation.cc)

install(TARGETS pk_prescribed_volumetric_deformation DESTINATION lib)

if ( BUILD_TESTS )

endif()
