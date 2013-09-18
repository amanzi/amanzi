# -*- mode: cmake -*-
#
# Define a project name
# After this command the following varaibles are defined
#   STRUCTUREDTRANS_SOURCE_DIR
#   STRUCTUREDTRANS_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(STRUCTUREDTRANS)

include_directories(${STRUCTUREDTRANS_SOURCE_DIR})

if (AMANZI_CHEMEVOL_PKG STREQUAL "AMANZI")
  include_directories(${CHEM_SOURCE_DIR}/includes)
  include_directories(${DBC_SOURCE_DIR})
endif()
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CCSE_INCLUDE_DIRS})
include_directories(${PMAMR_INCLUDE_DIRS})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})

set(CXX_header_files TensorDiffusionBndry.H TensorDiffusion_PK.H TensorOp_F.H TensorOp.H Diffuser.H MFVector.H ABecHelper.H)
set(CXX_source_files TensorDiffusionBndry.cpp TensorDiffusion_PK.cpp TensorOp.cpp MFVector.cpp ABecHelper.cpp)
set(F90_source_files )
set(F77_source_files )
set(FPP_source_files TensorOp_${AMANZI_SPACEDIM}D.F)

preprocess_boxlib_fortran(FPP_out_files ${FPP_source_files})

file(GLOB h_inc_files "*.H")
add_amanzi_library(strans
                   SOURCE ${FPP_out_files} ${CXX_source_files} ${F90_source_files} ${F77_source_files}
                   HEADERS ${h_inc_files})

add_amanzi_library(structtrans SOURCE ${CXX_source_files} HEADERS ${CXX_headers_files})

if (BUILD_TESTS)

  if (${CMAKE_SYSTEM_NAME} STREQUAL Darwin)
    #message(">>>>>> using darwin ld flags")
    if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
      set(ld_flag_prefix -force_load)
      set(ld_flag_suffix)
    endif()
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

  set(reqd_libs box_c box_cbndry box_cmg box_ctmg pmamr_sub strans)
  set(structured_grid_libs)
  foreach(lib ${reqd_libs})
    list(APPEND structured_grid_libs ${ld_flag_prefix} ${lib} ${ld_flag_suffix})
  endforeach(lib)
  # This is a hack to get a couple of the things in libbox_camr.a explicitly without including the whole thing...ick
  list(APPEND structured_grid_libs box_camr)

  if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU) 
    if (${CMAKE_SYSTEM_NAME} STREQUAL Linux)
      list(APPEND structured_grid_libs gfortran)
    endif()
  endif()

  set(TRANSTEST "${STRUCTUREDTRANS_SOURCE_DIR}/test")
  include_directories(${TRANSTEST})
  set(tpl_libs ${Teuchos_LIBRARIES})

  set(CXX_header_files_1 ${TRANSTEST}/tTensorDiffusion_PK_F.H)
  set(CXX_source_files_1 ${TRANSTEST}/tTensorDiffusion_PK.cpp)
  set(F90_source_files_1 )
  set(F77_source_files_1 )
  set(FPP_source_files_1 ${TRANSTEST}/tTensorDiffusion_PK_${AMANZI_SPACEDIM}D.F)

  preprocess_boxlib_fortran(FPP_out_files_1 ${FPP_source_files_1})

  add_amanzi_test(tTensorDiffusion_PK tTensorDiffusion_PK
                  KIND unit
                  SOURCE ${FPP_out_files_1} ${CXX_source_files_1} ${F90_source_files_1} ${F77_source_files_1}
                  LINK_LIBS ${structured_grid_libs} ${tpl_libs})

  if (NOT (CMAKE_SOURCE_DIR EQUAL CMAKE_BINARY_DIR))
    execute_process(COMMAND ${CMAKE_COMMAND} -E
                    copy ${TRANSTEST}/tTensorDiffusion_TestRes_H ${STRUCTUREDTRANS_BINARY_DIR})
    execute_process(COMMAND ${CMAKE_COMMAND} -E
                    copy ${TRANSTEST}/tTensorDiffusion_TestRes_D ${STRUCTUREDTRANS_BINARY_DIR})
  endif()
endif()

