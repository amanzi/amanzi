# -*- mode: cmake -*-
#
# Define a project name
# After this command the following varaibles are defined
#   SCHEM_SOURCE_DIR
#   SCHEM_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(SCHEM)

# Get the includes directory for chemistry
get_property(CHEM_INCLUDES_DIR GLOBAL PROPERTY CHEM_INCLUDES_DIR)

include_directories(${PMAMR_INCLUDE_DIRS})
include_directories(${ATK_SOURCE_DIR})
include_directories(${DBC_SOURCE_DIR})
include_directories(${PLATFORM_INTERFACE_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CCSE_INCLUDE_DIRS})
include_directories(${GEOCHEM_STATE_DIR}) 
include_directories(${GEOCHEM_SOURCE_DIR}) 
include_directories(${CHEM_INCLUDES_DIR})
include_directories(${ALQUIMIA_INCLUDE_DIRS})
add_definitions("-DENABLE_Structured")

# External (TPL) include directories
include_directories(${Trilinos_INCLUDE_DIRS})
include_directories(${PETSc_INCLUDE_DIRS})

set(CXX_header_files AlquimiaHelper_Structured.H AmanziChemHelper_Structured.H)

set(CXX_source_files ChemistryHelper_Structured.cpp AlquimiaHelper_Structured.cpp AmanziChemHelper_Structured.cpp)
set(F90_source_files )
set(F77_source_files )
set(FPP_source_files )

preprocess_boxlib_fortran(FPP_out_files ${FPP_source_files})

add_library(schem OBJECT ${FPP_out_files}
                         ${CXX_source_files}
                         ${CXX_header_files}
                         ${F90_source_files}
                         ${F77_source_files}
                         ${h_inc_files})

if (BUILD_TESTS)

  include_directories(${CCSE_INCLUDE_DIRS})
  link_directories(${CCSE_LIBRARY_DIR})

  set(structured_grid_libs structured_amanzi ${CCSE_LIBRARIES} ${PETSc_LIBRARIES})
  if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU) 
    if (${CMAKE_SYSTEM_NAME} STREQUAL Linux)
      list(APPEND structured_grid_libs gfortran)
    endif()
  endif()
  set(tpl_libs)

  get_property(GEOCHEM_LIBS GLOBAL PROPERTY GEOCHEM_LIBS)
  set(common_libs error_handling ${GEOCHEM_LIBS} interface atk)

  set(ssc_libs ${ssc_structured_grid_libs}
               ${ssc_unstructured_grid_libs}
               ${ssc_common_libs})
  set(tpl_libs ${amanzi_mesh_libs}
               ${Teuchos_LIBRARIES}
               ${Epetra_LIBRARIES}
               ${Boost_LIBRARIES}
               ${XERCES_LIBRARIES}
               ${XERCES_ICU_LIBRARIES})

  set(SCHEMTEST "${SCHEM_SOURCE_DIR}/test")

  #  tests removed #

endif()

