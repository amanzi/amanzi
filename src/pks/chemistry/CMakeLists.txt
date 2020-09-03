# -*- mode: cmake -*-

#
#  Amanzi
#    Chemistry Process Kernel
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)

#
# Define a project name
# After this command the following varaibles are defined
#   CHEMPK_SOURCE_DIR
#   CHEMPK_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(CHEMPK)

# Chemistry stuff.
get_property(CHEM_INCLUDES_DIR GLOBAL PROPERTY CHEM_INCLUDES_DIR)
get_property(GEOCHEM_LIBS GLOBAL PROPERTY GEOCHEM_LIBS)
include_directories(${CHEM_INCLUDES_DIR})

# External (TPL) include directories
# Although it does not call HDF5 directly 
# it does need HDF5 for the Vis.hpp file 
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${HDF5_INCLUDE_DIRS})
# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})
include_directories(${PKS_SOURCE_DIR})
include_directories(${GEOCHEM_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${MFUNCS_SOURCE_DIR})
include_directories(${MESH_DATA_SOURCE_DIR})
include_directories(${DATA_STRUCTURES_SOURCE_DIR})
include_directories(${FUNCS_SOURCE_DIR})
include_directories(${STATE_SOURCE_DIR})
include_directories(${OUTPUT_SOURCE_DIR})

#
# Install Targets
#
set(chem_inc_files 
                   Alquimia_PK.hh
                   Amanzi_PK.hh
                   Chemistry_PK.hh
                   pks_chemistry_registration.hh
)

#
# Library: chemistry_pk
#
set(chem_src_files)
if (ENABLE_Unstructured)
  list(APPEND chem_src_files Chemistry_PK.cc
                             Amanzi_PK.cc)
  if (ENABLE_ALQUIMIA)
    list(APPEND chem_src_files Alquimia_PK.cc)
  endif()
endif()

add_install_include_file(${chem_inc_files})

add_amanzi_library(chemistry_pk
                   SOURCE ${chem_src_files}
                   HEADERS ${chem_inc_files}
                   LINK_LIBS ${GEOCHEM_LIBS} pks ${Epetra_LIBRARIES}
                             ${XERCES_LIBRARIES} ${XERCES_ICU_LIBRARIES})

if (BUILD_TESTS AND ENABLE_Unstructured)
    # Add UnitTest include directory
    include_directories(${UnitTest_INCLUDE_DIRS})
    include_directories(${MESH_FACTORY_SOURCE_DIR})

    # Add the chemistry directory to the include paths
    include_directories(${CHEMPK_SOURCE_DIR})

    # Copy test directory files if any out of source build
    if (NOT (${CMAKE_SOURCE_DIR} EQUAL ${CMAKE_BINARY_DIR}))
        file(GLOB DataFiles "${CHEMPK_SOURCE_DIR}/test/*.xml"
                            "${CHEMPK_SOURCE_DIR}/test/*.bgd"
                            "${CHEMPK_SOURCE_DIR}/test/*.dbs"
                            "${CHEMPK_SOURCE_DIR}/test/*.in"
                            "${CHEMPK_SOURCE_DIR}/test/*.dat")
        file(COPY ${DataFiles} DESTINATION ${CHEMPK_BINARY_DIR}/test/)
    endif()

    add_amanzi_test(chemistry_amanzi_pk chemistry_amanzi_pk
                    KIND int
                    SOURCE test/Main.cc test/chemistry_amanzi_pk.cc
                    LINK_LIBS chemistry_pk state pks data_structures
                              mesh_functions mesh mesh_factory error_handling
                              ${UnitTest_LIBRARIES} ${Teuchos_LIBRARIES} ${Epetra_LIBRARIES}
                              ${XERCES_LIBRARIES} ${XERCES_ICU_LIBRARIES} ${NOX_LIBRARIES}
                              ${HYPRE_LIBRARY})

    if (ENABLE_ALQUIMIA)
      add_amanzi_test(chemistry_alquimia_pk chemistry_alquimia_pk
                      KIND int
                      SOURCE test/Main.cc test/chemistry_alquimia_pk.cc
                      LINK_LIBS chemistry_pk state pks data_structures
                                mesh_functions mesh mesh_factory error_handling
                                ${UnitTest_LIBRARIES} ${Teuchos_LIBRARIES} ${Epetra_LIBRARIES}
                                ${HDF5_LIBRARIES} ${XERCESC_LIBRARIES} ${XERCES_ICU_LIBRARIES}
                                ${NOX_LIBRARIES})
    endif()

    if (ENABLE_AmanziPhysicsModule)
      # this test uses the Amanzi executable...
      
      # The following test is disabled in parallel environments. It has to be rewritten in order 
      # to function properly there.
      add_amanzi_test(chemistry_benchmarks_1d chemistry_benchmark_1d
                      KIND int
                      SOURCE test/chemistry_benchmarks_1d.cc
                      LINK_LIBS chemistry_pk state pks data_structures
                              mesh_functions mesh mesh_factory error_handling
                              ${UnitTest_LIBRARIES} ${Teuchos_LIBRARIES} ${Epetra_LIBRARIES}
                              ${HDF5_LIBRARIES} ${XERCESC_LIBRARIES} ${XERCES_ICU_LIBRARIES}
                              ${NOX_LIBRARIES})
    endif()
endif()
