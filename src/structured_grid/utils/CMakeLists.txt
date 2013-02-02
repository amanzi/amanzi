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

set(CXX_header_files ParmParseHelpers.H InputParser_Structured.H
TabularFunction.H EventCoord.H WritePlotfile.H)
set(CXX_source_files ParmParseHelpers.cpp InputParser_Structured.cpp
TabularFunction.cpp EventCoord.cpp WritePlotfile.cpp)

add_amanzi_library(structutils SOURCE ${CXX_source_files} HEADERS ${CXX_headers_files})

if (BUILD_TESTS)

  if (${CMAKE_SYSTEM_NAME} STREQUAL Darwin)
    #message(">>>>>> using darwin ld flags")
    set(ld_flag_prefix -force_load)
    set(ld_flag_suffix)
  elseif (${CMAKE_SYSTEM_NAME} STREQUAL Linux)
    #message(">>>>>> using linux ld flags")
    set(ld_flag_prefix -Wl,-whole-archive)
    set(ld_flag_suffix -Wl,-no-whole-archive)
  else()
    message(FATAL_ERROR "Linking of structured grid is not supported on this system: '${CMAKE_SYSTEM_NAME}'")
  endif()

  include_directories(${CCSE_INCLUDE_DIRS})
  link_directories(${CCSE_LIBRARY_DIR})
  add_definitions("-DENABLE_Structured")

  set(reqd_libs structutils box_c)
  set(structured_grid_libs)
  foreach(lib ${reqd_libs})
    list(APPEND structured_grid_libs ${ld_flag_prefix} ${lib} ${ld_flag_suffix})
  endforeach(lib)

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

