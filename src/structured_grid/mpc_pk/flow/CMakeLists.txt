# -*- mode: cmake -*-
#
# Define a project name
# After this command the following varaibles are defined
#   STRUCTUREDFLOW_SOURCE_DIR
#   STRUCTUREDFLOW_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(STRUCTUREDFLOW)

include_directories(${STRUCTUREDFLOW_SOURCE_DIR})

if (AMANZI_CHEMEVOL_PKG STREQUAL "AMANZI")
  include_directories(${CHEM_SOURCE_DIR}/includes)
  include_directories(${DBC_SOURCE_DIR})
endif()
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

include_directories(${CHEM_SOURCE_DIR}/includes)
include_directories(${DBC_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CCSE_INCLUDE_DIRS})
include_directories(${STRUCTUREDFLOW_SOURCE_DIR})
include_directories(${STRUCTUREDUTILS_SOURCE_DIR})
include_directories(${PMAMR_INCLUDE_DIRS})
include_directories(${PETSC_INCLUDE_DIRS})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})

set(CXX_header_files RichardSolver.H RICHARDSOLVER_F.H NLScontrol.H NLSstatus.H RSdata.H RSAMRdata.H)
set(CXX_source_files RichardSolver.cpp NLScontrol.cpp RSdata.cpp RSAMRdata.cpp)
set(F90_source_files )
set(F77_source_files )
set(FPP_source_files RICHARDSOLVER_${AMANZI_SPACEDIM}D.F)

preprocess_boxlib_fortran(FPP_out_files ${FPP_source_files})

add_library(sflow OBJECT ${FPP_out_files}
                         ${CXX_source_files}
                         ${F90_source_files}
                         ${F77_source_files}
                         ${h_inc_files})

if (BUILD_TESTS)

  link_directories(${CCSE_LIBRARY_DIR})
  add_definitions("-DENABLE_Structured")

  if ( ENABLE_Unstructured )
    set(unstructured_grid_libs udriv mpc state transport flow chemistry_pk whetstone matrix_audit)
  else()
    set(unstructured_grid_libs)
  endif()

  # Note the following is rather icky...because of circular dependencies, these things have to be listed twice
  # This is similar to the main amanzi exec build, but worse.  Well....it works, for now.  Not sure what a better
  # long-term strategy is
  if(ENABLE_Structured)
    set(structured_grid_libs structured_amanzi ${CCSE_LIBRARIES} structured_amanzi ${CCSE_LIBRARIES} ${PETSC_LIBRARIES})
    if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU) 
      if (${CMAKE_SYSTEM_NAME} STREQUAL Linux)
        list(APPEND structured_grid_libs gfortran)
      endif()
    endif()
  else()
    set(structured_grid_libs)
  endif()

  get_property(GEOCHEM_LIBS GLOBAL PROPERTY GEOCHEM_LIBS)
  set(my_common_libs error_handling ${GEOCHEM_LIBS} interface atk)

  set(my_libs ${structured_grid_libs}
              ${unstructured_grid_libs}
              ${my_common_libs})
  set(my_tpl_libs ${amanzi_mesh_libs}
                  ${Teuchos_LIBRARIES}
                  ${Epetra_LIBRARIES}
                  ${Boost_LIBRARIES}
                  ${XERCES_LIBRARIES}
                  ${XERCES_ICU_LIBRARIES})

  if (ENABLE_Structured)
    list(APPEND my_tpl_libs ${CCSE_EXT_LIBRARIES})
  endif()


  set(FLOWTEST "${STRUCTUREDFLOW_SOURCE_DIR}/test")
  include_directories(${FLOWTEST})
  set(CXX_header_files_1 ${FLOWTEST}/RStdata.H)
  set(CXX_source_files_1 ${FLOWTEST}/tRS.cpp ${FLOWTEST}/RStdata.cpp)
  set(F90_source_files_1 )
  set(F77_source_files_1 )
  set(FPP_source_files_1 )

  preprocess_boxlib_fortran(FPP_out_files_1 ${FPP_source_files_1})

  add_amanzi_test(tRS tRS
                  KIND unit
                  SOURCE ${FPP_out_files_1} ${CXX_source_files_1} ${F90_source_files_1} ${F77_source_files_1}
                  LINK_LIBS ${my_libs} ${my_tpl_libs})

  if (NOT (CMAKE_SOURCE_DIR EQUAL CMAKE_BINARY_DIR))
    execute_process(COMMAND ${CMAKE_COMMAND} -E
                    copy ${FLOWTEST}/.petsc ${STRUCTUREDFLOW_BINARY_DIR})
  endif()
endif()
