# -*- mode: cmake -*-
#
# Define a project name
# After this command the following varaibles are defined
#   SADV_SOURCE_DIR
#   SADV_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(SADV)

include_directories(${SADV_SOURCE_DIR})

if (AMANZI_CHEMEVOL_PKG STREQUAL "AMANZI")
  include_directories(${DBC_SOURCE_DIR})
endif()
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CCSE_INCLUDE_DIRS})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})

include_directories(${PMAMR_INCLUDE_DIRS})

set(CXX_header_files Advection.H Advection_F.H)
set(CXX_source_files Advection.cpp)
set(F90_source_files Advection_${AMANZI_SPACEDIM}D.f90)
set(F77_source_files )
set(FPP_source_files )

preprocess_boxlib_fortran(FPP_out_files ${FPP_source_files})

add_library(sadv OBJECT ${FPP_out_files}
                        ${CXX_source_files}
                        ${F90_source_files}
                        ${F77_source_files})

