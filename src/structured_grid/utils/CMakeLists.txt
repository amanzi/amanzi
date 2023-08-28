# -*- mode: cmake -*-
#
# Define a project name
# After this command the following varaibles are defined
#   SUTILS_SOURCE_DIR
#   SUTILS_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(SUTILS)

include_directories(${SUTILS_SOURCE_DIR})

# we need to close all DSO objects
include_directories(${ATK_SOURCE_DIR})

if (AMANZI_CHEMEVOL_PKG STREQUAL "AMANZI")
  include_directories(${DBC_SOURCE_DIR})
endif()
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CCSE_INCLUDE_DIRS})
include_directories(${SAMR_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Trilinos_INCLUDE_DIRS})
include_directories(${PETSC_INCLUDE_DIRS})
include_directories(${PMAMR_INCLUDE_DIRS})

set(CXX_header_files ParmParseHelpers.H
                     InputParser_Structured.H
                     TabularFunction.H
                     EventCoord.H
                     WritePlotfile.H
                     MonotCubicInterpolator.H
                     MFTower_F.H
                     MFTower.H
                     Layout.H)

set(CXX_source_files ParmParseHelpers.cpp
                     InputParser_Structured.cpp
                     TabularFunction.cpp
                     EventCoord.cpp
                     WritePlotfile.cpp
                     MonotCubicInterpolator.cpp
                     MFTower.cpp
                     Layout.cpp)

set(FPP_source_files MFTower_${AMANZI_SPACEDIM}D.F)

preprocess_boxlib_fortran(FPP_out_files ${FPP_source_files})

add_library(sutils OBJECT ${FPP_out_files}
                        ${CXX_source_files}
                        ${F90_source_files}
                        ${F77_source_files}
                        ${CXX_header_files})

# This library depends on samr.
add_dependencies(sutils samr)

if (BUILD_TESTS)
  include_directories(${CCSE_INCLUDE_DIRS})
  link_directories(${CCSE_LIBRARY_DIR})

  set(structured_grid_libs structured_amanzi interface interface_s atk error_handling)

  if (ENABLE_ALQUIMIA)
    list(APPEND structured_grid_libs geochemistry)
  endif()

  if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU) 
    if (${CMAKE_SYSTEM_NAME} STREQUAL Linux)
      list(APPEND structured_grid_libs gfortran)
    endif()
  endif()

  set(tpl_libs ${Teuchos_LIBRARIES} ${Boost_LIBRARIES} 
               ${XERCES_LIBRARIES} ${XERCES_ICU_LIBRARIES}
               ${CCSE_LIBRARIES} ${PETSC_LIBRARIES} ${HDF5_LIBRARIES})

  if (BUILD_SHARED_LIBS)
    set(closing_dso ../shared_dso.f90)
  endif()

  add_amanzi_test(tPPH tPPH
                  KIND unit
                  SOURCE tParmParseHelpers.cpp ${closing_dso}
                  LINK_LIBS ${structured_grid_libs} ${tpl_libs})

  if (NOT (CMAKE_SOURCE_DIR EQUAL CMAKE_BINARY_DIR))
    execute_process(COMMAND ${CMAKE_COMMAND} -E
                    copy ${SUTILS_SOURCE_DIR}/sample.xml ${SUTILS_BINARY_DIR})
  endif()

  add_amanzi_test(tEC tEC
                  KIND unit
                  SOURCE tEventCoord.cpp ${closing_dso}
                  LINK_LIBS ${structured_grid_libs} ${tpl_libs})

endif()

