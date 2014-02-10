# -*- mode: cmake -*-------------------------------------------
include(TestManager)
include(PrintVariable)

#
# Define a project name
# After this command the following varaibles are defined
#   SREGION_SOURCE_DIR
#   SREGION_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(SREGION)


include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CCSE_INCLUDE_DIRS})
include_directories(${SREGION_INCLUDE_DIRS})
include_directories(${STRUCTUREDUTILS_SOURCE_DIR})

# External (TPL) include directories
#include_directories(${Teuchos_INCLUDE_DIRS})

set(CXX_source_files Region.cpp 
                     RegionManager.cpp 
                     RegionData.cpp)

set(CXX_header_files RegionData.H
                     Region.H 
                     RegionManager.H)

set(F90_source_files )
set(F77_source_files )
set(FPP_source_files )

preprocess_boxlib_fortran(FPP_out_files ${FPP_source_files})

add_library(sregion OBJECT ${FPP_out_files}
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

  #set(common_libs error_handling ${GEOCHEM_LIBS} interface atk)

  #set(tpl_libs ${Teuchos_LIBRARIES})

  add_amanzi_test(tRegion tRegion
                  KIND unit
                  SOURCE test/tRegion.cpp
                  LINK_LIBS ${structured_grid_libs} ${common_libs} ${tpl_libs})

endif()

