# -*- mode: cmake -*-
#
# Define a project name
# After this command the following varaibles are defined
#   SFLOW_SOURCE_DIR
#   SFLOW_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(SFLOW)

include_directories(${SFLOW_SOURCE_DIR})

if (AMANZI_CHEMEVOL_PKG STREQUAL "AMANZI")
  include_directories(${DBC_SOURCE_DIR})
endif()
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

include_directories(${DBC_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CCSE_INCLUDE_DIRS})
include_directories(${SFLOW_SOURCE_DIR})
include_directories(${SUTILS_SOURCE_DIR})  
include_directories(${SFLOW_SOURCE_DIR})  
include_directories(${STRANS_SOURCE_DIR})  
include_directories(${SPROB_SOURCE_DIR})  
include_directories(${SAMR_SOURCE_DIR})  
include_directories(${SCHEM_SOURCE_DIR})  
include_directories(${SREGION_SOURCE_DIR})   
include_directories(${SMATERIALS_SOURCE_DIR})   
include_directories(${PETSc_INCLUDE_DIRS})
include_directories(${CHEMPK_SOURCE_DIR})
include_directories(${GEOCHEM_SOURCE_DIR})
include_directories(${Alquimia_INCLUDE_DIRS})

# External (TPL) include directories
include_directories(${Trilinos_INCLUDE_DIRS})

include_directories(${PMAMR_INCLUDE_DIRS})

set(CXX_header_files RichardSolver.H RichardSolver_F.H NLScontrol.H NLSstatus.H RSdata.H RSAMRdata.H)
set(CXX_source_files RichardSolver.cpp NLScontrol.cpp RSdata.cpp RSAMRdata.cpp)
set(F90_source_files )
set(F77_source_files )
set(FPP_source_files RichardSolver_${AMANZI_SPACEDIM}D.F)

preprocess_boxlib_fortran(FPP_out_files ${FPP_source_files})

add_library(sflow OBJECT ${FPP_out_files}
                         ${CXX_source_files}
                         ${F90_source_files}
                         ${F77_source_files}
                         ${h_inc_files})

if (BUILD_TESTS)
  link_directories(${CCSE_LIBRARY_DIR})
  if (ENABLE_Unstructured)
    set(unstructured_grid_libs udriver mpc state transport flow whetstone)
  else()
    set(unstructured_grid_libs)
  endif()

  # Note the following is rather icky...because of circular dependencies, these things have to be listed twice
  # This is similar to the main amanzi exec build, but worse.  Well....it works, for now.  Not sure what a better
  # long-term strategy is
  set(structured_grid_libs structured_amanzi ${CCSE_LIBRARIES} 
                           structured_amanzi ${CCSE_LIBRARIES}
                           error_handling interface interface_s atk)

  if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU) 
    if (${CMAKE_SYSTEM_NAME} STREQUAL Linux)
      list(APPEND structured_grid_libs gfortran)
    endif()
  endif()

  get_property(GEOCHEM_LIBS GLOBAL PROPERTY GEOCHEM_LIBS)

  set(tpl_libs ${Teuchos_LIBRARIES} ${Epetra_LIBRARIES}
               ${CCSE_EXT_LIBRARIES} ${PETSc_LIBRARIES}
               ${XERCES_LIBRARIES} ${XERCES_ICU_LIBRARIES}
               ${Boost_LIBRARIES})

  if (BUILD_SHARED_LIBS)
    set(closing_dso ../../shared_dso.f90)
  endif()

  set(FLOWTEST "${SFLOW_SOURCE_DIR}/test")
  include_directories(${FLOWTEST})
  set(CXX_header_files_1 ${FLOWTEST}/RStdata.H)
  set(CXX_source_files_1 ${FLOWTEST}/tRS.cpp ${FLOWTEST}/RStdata.cpp)
  set(F90_source_files_1)
  set(F77_source_files_1)
  set(FPP_source_files_1)

  preprocess_boxlib_fortran(FPP_out_files_1 ${FPP_source_files_1})

  set(CXX_header_files_2 ${FLOWTEST}/RStdata.H)
  set(CXX_source_files_2 ${FLOWTEST}/tRS_steady.cpp ${FLOWTEST}/RStdata.cpp)
  set(F90_source_files_2 )
  set(F77_source_files_2 )
  set(FPP_source_files_2 )

  preprocess_boxlib_fortran(FPP_out_files_2 ${FPP_source_files_2})

  if(${AMANZI_SPACEDIM} EQUAL 2)
    # (seems to require redundant lib listing to resolve everything....not sure why
    add_amanzi_test(tRS tRS "tRS.in"
                    KIND unit
                    SOURCE ${FPP_out_files_1} ${CXX_source_files_1} ${F90_source_files_1} 
                           ${F77_source_files_1} ${closing_dso}
                    LINK_LIBS ${structured_grid_libs} ${GEOCHEM_LIBS} ${tpl_libs})

    add_amanzi_test(tRS_steady tRS_steady "tRS_steady.in"
                    KIND unit
                    SOURCE ${FPP_out_files_2} ${CXX_source_files_2} ${F90_source_files_2} 
                           ${F77_source_files_2} ${closing_dso}
		    LINK_LIBS ${structured_grid_libs} ${GEOCHEM_LIBS} ${tpl_libs})
  else()
    # (seems to require redundant lib listing to resolve everything....not sure why
    add_amanzi_test(tRS tRS "tRS_3D.in"
                    KIND unit
                    SOURCE ${FPP_out_files_1} ${CXX_source_files_1} ${F90_source_files_1} 
                           ${F77_source_files_1}
                    LINK_LIBS ${structured_grid_libs} ${tpl_libs})
  endif()

  if (NOT (CMAKE_SOURCE_DIR EQUAL CMAKE_BINARY_DIR))
  endif()

  # Copy input files if an out of source build
  if (NOT (${SFLOW_SOURCE_DIR} EQUAL ${SFLOW_BINARY_DIR}) )

    execute_process(COMMAND ${CMAKE_COMMAND} -E
                    copy ${FLOWTEST}/.petsc ${SFLOW_BINARY_DIR})

    if(${AMANZI_SPACEDIM} EQUAL 2)
      execute_process(COMMAND ${CMAKE_COMMAND} -E 
                      copy ${FLOWTEST}/tRS.in ${SFLOW_BINARY_DIR}) 
      execute_process(COMMAND ${CMAKE_COMMAND} -E 
                      copy ${FLOWTEST}/tRS_steady.in ${SFLOW_BINARY_DIR}) 
    else()
      execute_process(COMMAND ${CMAKE_COMMAND} -E 
                      copy ${FLOWTEST}/tRS_3D.in ${SFLOW_BINARY_DIR}) 
    endif()
  endif()

endif()
