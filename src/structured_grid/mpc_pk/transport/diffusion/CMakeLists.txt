# -*- mode: cmake -*-
#
# Define a project name
# After this command the following varaibles are defined
#   STRANS_SOURCE_DIR
#   STRANS_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(STRANS)

include_directories(${STRANS_SOURCE_DIR})

if (AMANZI_CHEMEVOL_PKG STREQUAL "AMANZI")
  include_directories(${CHEM_SOURCE_DIR}/includes)
  include_directories(${DBC_SOURCE_DIR})
endif()
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CCSE_INCLUDE_DIRS})
include_directories(${SAMR_SOURCE_DIR}) 
include_directories(${SFLOW_SOURCE_DIR}) 
include_directories(${STRANS_SOURCE_DIR}) 

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})

include_directories(${PMAMR_INCLUDE_DIRS})

set(CXX_header_files TensorDiffusionBndry.H
                     TensorDiffusion_PK.H
                     TensorOp_F.H
                     TensorOp.H
                     Diffuser.H
                     MFVector.H
                     ABecHelper.H)
set(CXX_source_files TensorDiffusionBndry.cpp
                     TensorDiffusion_PK.cpp
                     TensorOp.cpp
                     MFVector.cpp
                     ABecHelper.cpp)
set(F90_source_files )
set(F77_source_files )
set(FPP_source_files TensorOp_${AMANZI_SPACEDIM}D.F)

preprocess_boxlib_fortran(FPP_out_files ${FPP_source_files})

add_library(strans OBJECT ${FPP_out_files}
                          ${CXX_source_files}
                          ${F90_source_files}
                          ${F77_source_files})

if (BUILD_TESTS)
  include_directories(${ATK_SOURCE_DIR})
  include_directories(${CCSE_INCLUDE_DIRS})
  link_directories(${CCSE_LIBRARY_DIR})

  set(structured_grid_libs structured_amanzi interface interface_s atk)

  if (ENABLE_ALQUIMIA)
    list(APPEND structured_grid_libs geochemistry)
  endif()

  if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU) 
    if (${CMAKE_SYSTEM_NAME} STREQUAL Linux)
      list(APPEND structured_grid_libs gfortran)
    endif()
  endif()

  set(tpl_libs ${CCSE_LIBRARIES} ${PETSC_LIBRARIES}
               ${XERCES_LIBRARIES} ${XERCES_ICU_LIBRARIES}
               ${Boost_LIBRARIES})

  if (BUILD_SHARED_LIBS)
    set(closing_dso ../../../shared_dso.f90)
  endif()

  set(TRANSTEST "${STRANS_SOURCE_DIR}/test")
  include_directories(${TRANSTEST})
  set(CXX_header_files_1 ${TRANSTEST}/tTensorDiffusion_PK_F.H)
  set(CXX_source_files_1 ${TRANSTEST}/tTensorDiffusion_PK.cpp)
  set(F90_source_files_1 )
  set(F77_source_files_1 )
  set(FPP_source_files_1 ${TRANSTEST}/tTensorDiffusion_PK_${AMANZI_SPACEDIM}D.F)

  preprocess_boxlib_fortran(FPP_out_files_1 ${FPP_source_files_1})

  if(${AMANZI_SPACEDIM} EQUAL 2)
    add_amanzi_test(tTensorDiffusion_PK tTensorDiffusion_PK
                    KIND unit
                    SOURCE ${FPP_out_files_1} ${CXX_source_files_1} ${F90_source_files_1} 
                           ${F77_source_files_1} ${closing_dso}
                    LINK_LIBS ${structured_grid_libs} ${tpl_libs})
  endif()

  if (NOT (CMAKE_SOURCE_DIR EQUAL CMAKE_BINARY_DIR))
    if(${AMANZI_SPACEDIM} EQUAL 2)
      execute_process(COMMAND ${CMAKE_COMMAND} -E
                      copy ${TRANSTEST}/tTensorDiffusion_TestRes_H ${STRANS_BINARY_DIR})
      execute_process(COMMAND ${CMAKE_COMMAND} -E
                      copy ${TRANSTEST}/tTensorDiffusion_TestRes_D ${STRANS_BINARY_DIR})
    endif()
  endif()
endif()

