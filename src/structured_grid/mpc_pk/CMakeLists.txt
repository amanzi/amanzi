# -*- mode: cmake -*-------------------------------------------
include(TestManager)
include(PrintVariable)

#
# Define a project name
# After this command the following varaibles are defined
#   PMAMR_SOURCE_DIR
#   PMAMR_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(PMAMR)


include_directories(${CHEM_SOURCE_DIR}/includes)
include_directories(${DBC_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CCSE_INCLUDE_DIRS})
include_directories(${PMAMR_INCLUDE_DIRS})
include_directories(${PETSC_INCLUDE_DIRS})
include_directories(${CHEM_INCLUDES_DIR})
include_directories(${CHEMPK_SOURCE_DIR})

print_variable(CHEMPK_SOURCE_DIR)

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})

set(CXX_source_files Diffusion.cpp
                     Godunov.cpp
                     MacOperator.cpp
                     MacProj.cpp
                     Observation.cpp
                     PMBld.cpp
                     PM_setup.cpp
                     PorousMedia.cpp
                     Region.cpp
                     Source.cpp
                     ViscBndry.cpp
                     PressBndry.cpp
                     RegionData.cpp
                     PMAmr.cpp
                     Layout.cpp
                     MFTower.cpp)

set(CXX_header_files BC_TYPES_PM.H
                     DIFFUSION_F.H
                     Diffusion.H
                     GODCOMM_F.H
                     GODUNOV_F.H
                     Godunov.H
                     Layout.H
                     MACOPERATOR_F.H
                     MacOperator.H
                     MacOpMacDrivers.H
                     MACOUTFLOWBC_F.H
                     MACPROJ_F.H
                     MacProj.H
                     MFTOWER_F.H
                     MFTower.H
                     Observation.H
                     PMAmr.H
                     PMAMR_Labels.H
                     POROUS_F.H
                     POROUSMEDIA_F.H
                     PorousMedia.H
                     PressBndry.H
                     RegionData.H
                     Region.H
                     RegType.H
                     Source.H
                     ViscBndry.H
                     VISCOPERATOR_F.H)

set(F90_source_files )
set(F77_source_files )
set(FPP_source_files  EOS_F.F
                      GODUNOV_F.F
                      DIFFUSION_${AMANZI_SPACEDIM}D.F
                      GODUNOV_${AMANZI_SPACEDIM}D.F
                      MACOPERATOR_${AMANZI_SPACEDIM}D.F
                      MACPROJ_${AMANZI_SPACEDIM}D.F
                      POROUS_${AMANZI_SPACEDIM}D.F
                      POROUSMEDIA_${AMANZI_SPACEDIM}D.F
                      VISCOPERATOR_${AMANZI_SPACEDIM}D.F
                      MFTOWER_${AMANZI_SPACEDIM}D.F)

preprocess_boxlib_fortran(FPP_out_files ${FPP_source_files})

add_library(pmamr OBJECT ${FPP_out_files}
                         ${CXX_source_files}
                         ${F90_source_files}
                         ${F77_source_files}
                         ${CXX_header_files})

if (BUILD_TESTS)

  include_directories(${CCSE_INCLUDE_DIRS})
  link_directories(${CCSE_LIBRARY_DIR})
  add_definitions("-DENABLE_Structured")

  set(structured_grid_libs structured_amanzi ${PETSC_LIBRARIES} ${CCSE_LIBRARIES})
  if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
    if (${CMAKE_SYSTEM_NAME} STREQUAL Linux)
      list(APPEND structured_grid_libs gfortran)
    endif()
  endif()

  get_property(GEOCHEM_LIBS GLOBAL PROPERTY GEOCHEM_LIBS)
  set(common_libs error_handling ${GEOCHEM_LIBS} interface atk)

  set(tpl_libs ${Teuchos_LIBRARIES})

  add_amanzi_test(tRegion tRegion
                  KIND unit
                  SOURCE tRegion.cpp
                  LINK_LIBS ${structured_grid_libs} ${common_libs} ${tpl_libs})

endif()

