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


if (AMANZI_CHEMEVOL_PKG STREQUAL "AMANZI")
  include_directories(${CHEM_SOURCE_DIR}/includes)
  include_directories(${DBC_SOURCE_DIR})
endif()
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CCSE_INCLUDE_DIRS})
include_directories(${PMAMR_INCLUDE_DIRS})

set(CXX_source_files Diffusion.cpp Godunov.cpp MacOperator.cpp MacProj.cpp Observation.cpp PMBld.cpp PM_setup.cpp PorousMedia.cpp Region.cpp Rock.cpp Source.cpp ViscBndry.cpp PressBndry.cpp RegionData.cpp PMAmr.cpp Layout.cpp MFTower.cpp Richard.cpp)
set(F90_source_files bivar.f90)
set(F77_source_files math_d.f vode.f)

set(FPP_source_files DERIVE_${AMANZI_SPACEDIM}D.F DIFFUSION_${AMANZI_SPACEDIM}D.F EOS_${AMANZI_SPACEDIM}D.F EOS_F.F GODUNOV_${AMANZI_SPACEDIM}D.F GODUNOV_F.F MACOPERATOR_${AMANZI_SPACEDIM}D.F MACPROJ_${AMANZI_SPACEDIM}D.F POROUS_${AMANZI_SPACEDIM}D.F POROUSMEDIA_${AMANZI_SPACEDIM}D.F TURN_F.F VISCOPERATOR_${AMANZI_SPACEDIM}D.F MFTOWER_${AMANZI_SPACEDIM}D.F RICHARD_${AMANZI_SPACEDIM}D.F)

if (AMANZI_CHEMEVOL_PKG STREQUAL "COREREACT")
  list(APPEND CXX_source_files PChemDriver.cpp)
  list(APPEND FPP_source_files CHEMDRIVER_${AMANZI_SPACEDIM}D.F)
  list(APPEND F90_source_files ChemModules.f90 ChemAllocate.f90 ChemBatch.f90 ChemRead.f90 ChemVode.f90)
endif()

preprocess_boxlib_fortran(FPP_out_files ${FPP_source_files})

add_library(pmamr ${FPP_out_files} ${CXX_source_files} ${F90_source_files} ${F77_source_files})
add_install_library(pmamr)

file(GLOB _inc_files "*.H")
add_install_include_file(${_inc_files})

if (BUILD_TESTS)
endif()

