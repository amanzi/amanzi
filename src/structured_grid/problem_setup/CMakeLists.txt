# -*- mode: cmake -*-------------------------------------------
#
# Define a project name
# After this command the following varaibles are defined
#   SPROB_SOURCE_DIR
#   SPROB_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(SPROB)

if (AMANZI_CHEMEVOL_PKG STREQUAL "AMANZI")
  include_directories(${DBC_SOURCE_DIR})
endif()
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CCSE_INCLUDE_DIRS})
include_directories(${SAMR_SOURCE_DIR})

include_directories(${PMAMR_INCLUDE_DIRS})

set(CXX_source_files)
set(F90_source_files)
set(F77_source_files)
set(FPP_source_files Prob_${AMANZI_SPACEDIM}D.F)

preprocess_boxlib_fortran(FPP_out_files ${FPP_source_files})

add_library(sprob OBJECT ${FPP_out_files} ${CXX_source_files} ${F90_source_files} ${F77_source_files})

if (BUILD_TESTS)
endif()

