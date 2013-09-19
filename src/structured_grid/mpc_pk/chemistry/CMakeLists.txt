# -*- mode: cmake -*-
#
# Define a project name
# After this command the following varaibles are defined
#   SCHEM_SOURCE_DIR
#   SCHEM_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(SCHEM)

include_directories(${SCHEM_SOURCE_DIR})

if (AMANZI_CHEMEVOL_PKG STREQUAL "AMANZI")
  include_directories(${CHEM_SOURCE_DIR}/includes)
  include_directories(${DBC_SOURCE_DIR})
endif()
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CCSE_INCLUDE_DIRS})

set(CXX_header_files MultiFabCellIterator.H
                     MultiFabCellVector.H)

set(CXX_source_files MultiFabCellIterator.cpp)
set(F90_source_files )
set(F77_source_files )
set(FPP_source_files )

preprocess_boxlib_fortran(FPP_out_files ${FPP_source_files})

add_library(schem OBJECT ${FPP_out_files}
                         ${CXX_source_files}
                         ${F90_source_files}
                         ${F77_source_files}
                         ${h_inc_files})

if (BUILD_TESTS)

  include_directories(${CCSE_INCLUDE_DIRS})
  link_directories(${CCSE_LIBRARY_DIR})
  add_definitions("-DENABLE_Structured")

  set(structured_grid_libs structured_amanzi ${CCSE_LIBRARIES})
  if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU) 
    if (${CMAKE_SYSTEM_NAME} STREQUAL Linux)
      list(APPEND structured_grid_libs gfortran)
    endif()
  endif()
  set(tpl_libs)

  set(SCHEMTEST "${SCHEM_SOURCE_DIR}/test")
  add_amanzi_test(tMultiFabCellVector tMultiFabCellVector
                  KIND unit
                  SOURCE ${SCHEMTEST}/tMultiFabCellVector.cpp
                  LINK_LIBS ${structured_grid_libs} ${tpl_libs})
endif()

