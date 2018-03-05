# -*- mode: cmake -*-
#
# Define a project name
# After this command the following varaibles are defined
#   SMATERIALS_SOURCE_DIR
#   SMATERIALS_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(SMATERIALS)

set(CXX_header_files Property.H
                     Material.H
                     MatFiller.H
                     MatFiller_F.H
                     RockManager.H
                     RockManager_F.H
                     GSLibInt.H
                     GSLibInt_F.H)

set(CXX_source_files Property.cpp
                     Material.cpp
                     MatFiller.cpp
                     RockManager.cpp
                     GSLibInt.cpp)

set(F90_source_files sgsim.f90
                     geostat.f90)

if(${AMANZI_SPACEDIM} EQUAL 3)
  list (APPEND F90_source_files bivar.f90)
endif()

set(FPP_source_files MatFiller_${AMANZI_SPACEDIM}D.F
                     GSLibInt_${AMANZI_SPACEDIM}D.F
                     RockManager_${AMANZI_SPACEDIM}D.F)

include_directories(${SREGION_SOURCE_DIR})  
include_directories(${SMATERIALS_SOURCE_DIR})  
include_directories(${SAMR_SOURCE_DIR})  
include_directories(${ATK_SOURCE_DIR})
include_directories(${PLATFORM_INTERFACE_SOURCE_DIR})

# The interfaces
include_directories(${INTERFACE_PLATFORM_DIR})
include_directories(${AMANZI_BINARY_DIR}) # required to pick up amanzi_version.hh

#
# The structured grid pieces
#
include_directories(${CCSE_INCLUDE_DIRS})
include_directories(${SAMR_INCLUDE_DIRS})
link_directories(${CCSE_LIBRARY_DIR})
include_directories(${SDRIVER_SOURCE_DIR})
link_directories(${SDRIVER_SOURCE_DIR})
include_directories(${SUTILS_SOURCE_DIR})
link_directories(${SUTILS_SOURCE_DIR})
add_definitions("-DENABLE_Structured")

#
# The common pieces
#
include_directories(${DBC_SOURCE_DIR})
include_directories(${CHEM_INCLUDES_DIR})

# External (TPL) include directories
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${PETSC_INCLUDE_DIRS})
include_directories(${PMAMR_INCLUDE_DIRS})

preprocess_boxlib_fortran(FPP_out_files ${FPP_source_files})

add_library(smaterials OBJECT ${FPP_out_files}
                              ${CXX_source_files}
                              ${F90_source_files}
                              ${F77_source_files}
                              ${CXX_headers_files})

if (BUILD_TESTS) 
  include_directories(${CCSE_INCLUDE_DIRS})
  link_directories(${CCSE_LIBRARY_DIR})
  add_definitions("-DENABLE_Structured")

  set(structured_grid_libs structured_amanzi interface interface_s atk error_handling)

  if (ENABLE_ALQUIMIA)
    list(APPEND structured_grid_libs geochemistry)
  endif()

  if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
    if (${CMAKE_SYSTEM_NAME} STREQUAL Linux)
      list(APPEND structured_grid_libs gfortran)
    endif()
  endif()

  set(tpl_libs ${CCSE_LIBRARIES} ${PETSC_LIBRARIES} ${CCSE_LIBRARIES}
               ${XERCES_LIBRARIES} ${XERCES_ICU_LIBRARIES}
               ${Boost_LIBRARIES} ${HDF5_LIBRARIES})

  include_directories(${SMATERIALS_SOURCE_DIR}/test) 

  if (BUILD_SHARED_LIBS)
    set(closing_dso ../shared_dso.f90)
  endif()

  add_amanzi_test(RockManager tRM "tRM.in"
                  KIND unit
                  SOURCE test/tRM.cpp ${closing_dso}
                  LINK_LIBS ${structured_grid_libs} ${tpl_libs})

  add_amanzi_test(RockManager_fort tRM_fort "tRM.in"
                  KIND unit
                  SOURCE test/tRM_fort.cpp test/tRM_fort_${AMANZI_SPACEDIM}D.F ${closing_dso}
                  LINK_LIBS ${structured_grid_libs} ${tpl_libs})

  if(${AMANZI_SPACEDIM} EQUAL 2)
    add_amanzi_test(MaterialIDFiller tMIDf
                    KIND unit
                    SOURCE test/tMIDf.cpp ${closing_dso}
                    LINK_LIBS ${structured_grid_libs} ${tpl_libs})

    add_amanzi_test(MIDfill_RM tMIDf_RM "tMIDf_RM.in"
                    KIND unit
                    SOURCE test/tMIDf_RM.cpp ${closing_dso}
                    LINK_LIBS ${structured_grid_libs} ${tpl_libs})

    add_amanzi_test(RockManager_gs tRM_GS "tRM_GS.in" "nLevs=1"
                    KIND unit
                    SOURCE test/tRM_GS.cpp ${closing_dso}
                    LINK_LIBS ${structured_grid_libs} ${tpl_libs})

    add_amanzi_test(RockManager_gs2 RockManager_gs_parallel "./tRM_GS.in"
                    KIND unit
                    NPROCS 2
                    SOURCE test/tRM_GS.cpp ${closing_dso}
                    LINK_LIBS ${structured_grid_libs} ${tpl_libs})

    add_amanzi_test(MaterialPropertyFiller tMatFiller
                    KIND unit
                    SOURCE test/tMatFiller.cpp ${closing_dso}
                    LINK_LIBS ${structured_grid_libs} ${tpl_libs})

    add_amanzi_test(MaterialPropertyFillerC tMatFillerC
                    KIND unit
                    SOURCE test/tMatFillerC.cpp ${closing_dso}
                    LINK_LIBS ${structured_grid_libs} ${tpl_libs})

    add_amanzi_test(MaterialPropertyFiller4 tMatFiller_parallel "case_size=large"
                    KIND unit
                    NPROCS 4
                    SOURCE test/tMatFiller.cpp ${closing_dso}
                    LINK_LIBS ${structured_grid_libs} ${tpl_libs})
  endif()

  # Copy input files if an out of source build
  if (NOT (${SMATERIALS_SOURCE_DIR} EQUAL ${SMATERIALS_BINARY_DIR}) )

    execute_process(COMMAND ${CMAKE_COMMAND} -E 
      copy ${SMATERIALS_SOURCE_DIR}/test/tRM.in ${SMATERIALS_BINARY_DIR}) 

    if(${AMANZI_SPACEDIM} EQUAL 2)
      execute_process(COMMAND ${CMAKE_COMMAND} -E 
        copy ${SMATERIALS_SOURCE_DIR}/test/tMIDf_RM.in ${SMATERIALS_BINARY_DIR}) 
      execute_process(COMMAND ${CMAKE_COMMAND} -E 
        copy ${SMATERIALS_SOURCE_DIR}/test/tRM_GS.in ${SMATERIALS_BINARY_DIR}) 
      execute_process(COMMAND ${CMAKE_COMMAND} -E 
        copy ${SMATERIALS_SOURCE_DIR}/test/cluster.dat ${SMATERIALS_BINARY_DIR}) 
      execute_process(COMMAND ${CMAKE_COMMAND} -E 
        copy ${SMATERIALS_SOURCE_DIR}/test/sgsim.par ${SMATERIALS_BINARY_DIR}) 
    endif()
  endif()

endif()
