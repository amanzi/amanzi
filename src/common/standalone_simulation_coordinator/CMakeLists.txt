# -*- mode: cmake -*-
#
# Define a project name
# After this command the following varaibles are defined
#   SSC_SOURCE_DIR
#   SSC_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.

project(SSC)

if ( ENABLE_Structured)
  include(CCSEOptions)
endif()

include(BinaryManager)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})
include_directories(${PLATFORM_INTERFACE_SOURCE_DIR})

include_directories(${AMANZI_BINARY_DIR}) # required to pick up amanzi_version.hh

#
# The structured grid pieces
#
if ( ENABLE_Structured )
  include_directories(${CCSE_INCLUDE_DIRS})
  include_directories(${SAMR_INCLUDE_DIRS})
  link_directories(${CCSE_LIBRARY_DIR})
  include_directories(${SDRIVER_SOURCE_DIR})
  link_directories(${SDRIVER_SOURCE_DIR})
  include_directories(${SUTILS_SOURCE_DIR})
  link_directories(${SUTILS_SOURCE_DIR})
  add_definitions("-DENABLE_Structured")
endif()

#
# The common pieces
#
include_directories(${DBC_SOURCE_DIR})
include_directories(${CHEM_INCLUDES_DIR})

#
# The unstructured grid pieces
#
if ( ENABLE_Unstructured )
  include_directories(${ATK_SOURCE_DIR})
  include_directories(${FUNCS_SOURCE_DIR})
  include_directories(${MESH_FACTORY_SOURCE_DIR})
  include_directories(${MESH_SOURCE_DIR})
  include_directories(${MESH_DATA_SOURCE_DIR})
  include_directories(${MFUNCS_SOURCE_DIR})
  include_directories(${DATA_STRUCTURES_SOURCE_DIR})
  include_directories(${GEOMETRY_SOURCE_DIR})
  include_directories(${OPERATORS_SOURCE_DIR})
  include_directories(${OUTPUT_SOURCE_DIR})
  include_directories(${SOLVERS_SOURCE_DIR})
  include_directories(${STATE_SOURCE_DIR})
  include_directories(${STATE_SOURCE_DIR}/data)
  include_directories(${TRANSPORT_SOURCE_DIR})
  include_directories(${FLOW_SOURCE_DIR})
  include_directories(${MPC_SOURCE_DIR})
  include_directories(${EXO_SOURCE_DIR})
  include_directories(${EXEC_SOURCE_DIR})
  include_directories(${WHETSTONE_SOURCE_DIR})

  link_directories(${EXEC_SOURCE_DIR})

# The mesh frameworks
  set(amanzi_mesh_libs "mesh_factory")

  add_definitions("-DENABLE_Unstructured")
endif()

# External (TPL) include directories
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${XERCES_INCLUDE_DIRS})
include_directories(${Alquimia_INCLUDE_DIRS})
include_directories(${HDF5_INCLUDE_DIRS})

#
# More unstructured grid pieces
#
if (ENABLE_Unstructured)
endif()


#
# Binary: amanzi-demo-driver
#
if (ENABLE_Unstructured)
  set(ssc_unstructured_grid_libs udriv mpc state transport flow energy pks eos operators whetstone chemistry_pk)
else()
  set(ssc_unstructured_grid_libs)
endif()

if(ENABLE_Structured)
  set(ssc_structured_grid_libs interface_s structured_amanzi ${CCSE_LIBRARIES} structured_amanzi)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU) 
    if (${CMAKE_SYSTEM_NAME} STREQUAL Linux)
      list(APPEND ssc_structured_grid_libs gfortran)
    endif()
  endif()
else()
  set(ssc_structured_grid_libs)
endif()

get_property(GEOCHEM_LIBS GLOBAL PROPERTY GEOCHEM_LIBS)
set(ssc_common_libs error_handling ${GEOCHEM_LIBS} interface atk)

set(ssc_libs ${ssc_structured_grid_libs}
             ${ssc_unstructured_grid_libs}
             ${ssc_common_libs})
set(ssc_tpl_libs ${amanzi_mesh_libs}
                 ${Teuchos_LIBRARIES}
                 ${Epetra_LIBRARIES}
                 ${Boost_LIBRARIES}
                 ${XERCES_LIBRARIES}
                 ${XERCES_ICU_LIBRARIES})

if (ENABLE_Structured)
  list(APPEND ssc_tpl_libs ${CCSE_EXT_LIBRARIES})
      if (BUILD_SHARED_LIBS)
        set(closing_dso main_s_dso.f90)
      endif()
endif()

list(APPEND ssc_tpl_libs ${PETSc_LIBRARIES})

# 
# Main: Amanzi
#
add_amanzi_executable(amanzi
                      SOURCE Main.cc SimulatorFactory.cc ${closing_dso}
                      LINK_LIBS ${ssc_libs} ${ssc_tpl_libs} 
                      OUTPUT_NAME amanzi
                      OUTPUT_DIRECTORY ${SSC_BINARY_DIR})

if (BUILD_TESTS) 
    # Copy test directory files if an out of source build
    if (NOT (${SSC_SOURCE_DIR} EQUAL ${SSC_BINARY_DIR}) )
      execute_process(COMMAND ${CMAKE_COMMAND} -E 
          copy_directory ${SSC_SOURCE_DIR}/test ${SSC_BINARY_DIR}/test) 
    endif()

    if (ENABLE_Structured)
      # add_amanzi_test(ssc-structured-transport-0 amanzi "--xml_file=test/t0-structured-isv2.xml" KIND unit)
      add_amanzi_test(ssc-structured-flow-1 amanzi "--xml_file=test/t1-structured-isv2.xml" KIND unit)
    endif()

    if (ENABLE_Unstructured)
      add_amanzi_test(ssc_chemistry_amanzi amanzi "--xml_file=test/ssc_chemistry_amanzi.xml" KIND unit)
    endif()
endif()

