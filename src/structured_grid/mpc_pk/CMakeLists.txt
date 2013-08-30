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

set(CXX_source_files Diffusion.cpp Godunov.cpp MacOperator.cpp
MacProj.cpp Observation.cpp PMBld.cpp PM_setup.cpp PorousMedia.cpp
Region.cpp Source.cpp ViscBndry.cpp PressBndry.cpp
RegionData.cpp PMAmr.cpp Layout.cpp MFTower.cpp RichardSolver.cpp)
set(F90_source_files )
set(F77_source_files )

set(FPP_source_files 
DIFFUSION_${AMANZI_SPACEDIM}D.F EOS_F.F
GODUNOV_${AMANZI_SPACEDIM}D.F GODUNOV_F.F
MACOPERATOR_${AMANZI_SPACEDIM}D.F MACPROJ_${AMANZI_SPACEDIM}D.F
POROUS_${AMANZI_SPACEDIM}D.F POROUSMEDIA_${AMANZI_SPACEDIM}D.F
VISCOPERATOR_${AMANZI_SPACEDIM}D.F
MFTOWER_${AMANZI_SPACEDIM}D.F RICHARDSOLVER_${AMANZI_SPACEDIM}D.F)

preprocess_boxlib_fortran(FPP_out_files ${FPP_source_files})

file(GLOB h_inc_files "*.H")
add_amanzi_library(pmamr
                   SOURCE ${FPP_out_files} ${CXX_source_files} ${F90_source_files} ${F77_source_files}
                   HEADERS ${h_inc_files})

add_amanzi_library(pmamr_sub
                   SOURCE Region.cpp
                   HEADERS ${h_inc_files})


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

  set(reqd_libs pmamr_sub structutils box_c)
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

  add_amanzi_test(tRegion tRegion
                  KIND unit
                  SOURCE tRegion.cpp
                  LINK_LIBS ${structured_grid_libs} ${tpl_libs})

endif()

