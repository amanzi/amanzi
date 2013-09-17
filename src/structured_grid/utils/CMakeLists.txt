# -*- mode: cmake -*-
#
# Define a project name
# After this command the following varaibles are defined
#   STRUCTUREDUTILS_SOURCE_DIR
#   STRUCTUREDUTILS_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(STRUCTUREDUTILS)

include_directories(${STRUCTUREDUTILS_SOURCE_DIR})


if (AMANZI_CHEMEVOL_PKG STREQUAL "AMANZI")
  include_directories(${CHEM_SOURCE_DIR}/includes)
  include_directories(${DBC_SOURCE_DIR})
endif()
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CCSE_INCLUDE_DIRS})
include_directories(${PMAMR_INCLUDE_DIRS})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})

set(CXX_header_files ParmParseHelpers.H
                     InputParser_Structured.H
                     TabularFunction.H
                     EventCoord.H
                     WritePlotfile.H)
set(CXX_source_files ParmParseHelpers.cpp
                     InputParser_Structured.cpp
                     TabularFunction.cpp
                     EventCoord.cpp
                     WritePlotfile.cpp)

add_library(structutils OBJECT ${CXX_source_files}
                               ${CXX_headers_files})

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

  set(tpl_libs ${Teuchos_LIBRARIES})

  add_amanzi_test(tPPH tPPH
                  KIND unit
                  SOURCE tParmParseHelpers.cpp
                  LINK_LIBS ${structured_grid_libs} ${tpl_libs})

  if (NOT (CMAKE_SOURCE_DIR EQUAL CMAKE_BINARY_DIR))
    execute_process(COMMAND ${CMAKE_COMMAND} -E
                    copy ${STRUCTUREDUTILS_SOURCE_DIR}/sample.xml ${STRUCTUREDUTILS_BINARY_DIR})
  endif()

  add_amanzi_test(tEC tEC
                  KIND unit
                  SOURCE tEventCoord.cpp
                  LINK_LIBS ${structured_grid_libs} ${tpl_libs})

endif()

