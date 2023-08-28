# -*- mode: cmake -*-------------------------------------------
include(TestManager)
include(PrintVariable)

#
# Define a project name
# After this command the following varaibles are defined
#   SAMR_SOURCE_DIR
#   SAMR_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(SAMR)

# Get the includes directory for chemistry
get_property(CHEM_INCLUDES_DIR GLOBAL PROPERTY CHEM_INCLUDES_DIR)

include_directories(${DBC_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CCSE_INCLUDE_DIRS})
include_directories(${PMAMR_INCLUDE_DIRS})

include_directories(${PETSc_INCLUDE_DIRS})
include_directories(${CHEM_INCLUDES_DIR})
include_directories(${CHEMPK_SOURCE_DIR})
include_directories(${GEOCHEM_STATE_DIR}) 
include_directories(${GEOCHEM_SOURCE_DIR})
include_directories(${CHEM_INCLUDES_DIR})
include_directories(${Alquimia_INCLUDE_DIRS})

# External (TPL) include directories
include_directories(${Trilinos_INCLUDE_DIRS})
include_directories(${HDF5_INCLUDE_DIRS})

set(CXX_source_files Observation.cpp
                     PMBld.cpp
                     PM_setup.cpp
                     PorousMedia.cpp
                     Source.cpp
                     ViscBndry.cpp
                     PMAmr.cpp)

set(CXX_header_files BC_TYPES_PM.H
                     Observation.H
                     PMAmr.H
                     PMAMR_Labels.H
                     PorousMedia_F.H
                     PorousMedia.H
                     RegType.H
                     Source.H
                     ViscBndry.H
                     probdata.H)

set(F90_source_files )
set(F77_source_files )
set(FPP_source_files  PorousMedia_${AMANZI_SPACEDIM}D.F
                      Eos_F.F)

preprocess_boxlib_fortran(FPP_out_files ${FPP_source_files})

add_library(samr OBJECT ${FPP_out_files}
                        ${CXX_source_files}
                        ${F90_source_files}
                        ${F77_source_files}
                        ${CXX_header_files})

if (BUILD_TESTS)

  include_directories(${CCSE_INCLUDE_DIRS})
  link_directories(${CCSE_LIBRARY_DIR})

  set(structured_grid_libs structured_amanzi ${PETSc_LIBRARIES} ${CCSE_LIBRARIES})
  if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
    if (${CMAKE_SYSTEM_NAME} STREQUAL Linux)
      list(APPEND structured_grid_libs gfortran)
    endif()
  endif()

  get_property(GEOCHEM_LIBS GLOBAL PROPERTY GEOCHEM_LIBS)
  set(common_libs error_handling ${GEOCHEM_LIBS} interface atk)

  set(tpl_libs ${Teuchos_LIBRARIES} ${HDF5_LIBRARIES})

  #add_amanzi_test(tRegion tRegion
  #                KIND unit
  #                SOURCE tRegion.cpp
  #                LINK_LIBS ${structured_grid_libs} ${common_libs} ${tpl_libs})

endif()

