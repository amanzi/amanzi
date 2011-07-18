# -*- mode: cmake -*-
#

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CCSE_INCLUDE_DIRS})
include_directories(${PMAMR_INCLUDE_DIRS})

set(CXX_source_files PGslib.cpp)
set(F90_source_files geostat.f90 sgsim.f90)
set(F77_source_files)
set(FPP_files GSLIB_${AMANZI_SPACEDIM}D.F)

preprocess_boxlib_fortran(FPP_out_files ${FPP_files})

add_library(gsint ${FPP_out_files} ${CXX_source_files} ${F90_source_files} ${F77_source_files})
install_library(gsint)
link_list_add(gsint)

if (BUILD_TESTS)

endif()

