# -*- mode: cmake -*-

#
#  Amanzi
#    Amanzi unstructured driver
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)
include(LibraryManager)

#
# Define a project name
# After this command the following varaibles are defined
#   UDRIVER_SOURCE_DIR
#   UDRIVER_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(UDRIVER)

# Amanzi include directories
include_directories(${ATK_SOURCE_DIR})
include_directories(${CHEM_SOURCE_DIR}/base_chemistry)
include_directories(${CHEM_SOURCE_DIR}/reactions)
include_directories(${CHEM_SOURCE_DIR}/solvers)
include_directories(${CHEM_SOURCE_DIR}/utility)
include_directories(${CHEMPK_SOURCE_DIR})
include_directories(${DATA_STRUCTURES_SOURCE_DIR})
include_directories(${DBC_SOURCE_DIR})
include_directories(${DBG_SOURCE_DIR})
include_directories(${ENERGY_SOURCE_DIR})
include_directories(${EOS_SOURCE_DIR})
include_directories(${FLOW_SOURCE_DIR})
include_directories(${FUNCS_SOURCE_DIR})
include_directories(${GEOCHEM_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${MESH_EXTRACTED_SOURCE_DIR})
include_directories(${MESH_FACTORY_SOURCE_DIR})
include_directories(${MFUNCS_SOURCE_DIR})
include_directories(${MPC_SOURCE_DIR})
include_directories(${MPC_TREE_SOURCE_DIR})
include_directories(${MULTIPHASE_SOURCE_DIR})
include_directories(${OPERATORS_SOURCE_DIR})
include_directories(${OUTPUT_SOURCE_DIR})
include_directories(${PKS_SOURCE_DIR})
include_directories(${PLATFORM_INTERFACE_SOURCE_DIR})
include_directories(${SHALLOW_WATER_SOURCE_DIR})
include_directories(${SOLVERS_SOURCE_DIR})
include_directories(${SSC_SOURCE_DIR})
include_directories(${STATE_SOURCE_DIR})
include_directories(${STATE_SOURCE_DIR}/data)
include_directories(${STATE_SOURCE_DIR}/evaluators)
include_directories(${TIME_INTEGRATION_SOURCE_DIR})
include_directories(${TRANSPORT_SOURCE_DIR})
include_directories(${UDRIVER_SOURCE_DIR})
include_directories(${WHETSTONE_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
if(ENABLE_MESH_MOAB)
    include_directories(${MOAB_INCLUDE_DIRS})
    include_directories(${MOAB_SOURCE_DIR})
endif()
include_directories(${HDF5_C_INCLUDE_DIR})
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${ASCEMIO_INCLUDE_DIR})

# library
  add_amanzi_library(udriver
                   SOURCE AmanziUnstructuredGridSimulationDriver.cc
                   LINK_LIBS mesh mesh_factory state mpc eos 
                             flow transport energy shallow_water multiphase
                             interface error_handling
                             ${Epetra_LIBRARIES} ${Teuchos_LIBRARIES} ${Boost_LIBRARIES}
                             ${XERCES_LIBRARIES} ${XERCES_ICU_LIBRARIES})

add_install_include_file(AmanziUnstructuredGridSimulationDriver.hh)

if (BUILD_TESTS) 
    # Add UnitTest include directory
    include_directories(${UnitTest_INCLUDE_DIRS})

    # Copy test directory files if any out of source build
    if (NOT (${UDRIVER_SOURCE_DIR} EQUAL ${UDRIVER_BINARY_DIR}) )
        file(GLOB DataFiles "${UDRIVER_SOURCE_DIR}/test/*.xml"
                            "${UDRIVER_SOURCE_DIR}/test/*.exo")
        file(COPY ${DataFiles} DESTINATION ${UDRIVER_BINARY_DIR}/test/)
    endif()

    set(amanzi_libs geometry mesh mesh_factory state udriver flow pks operators whetstone data_structures)

    add_amanzi_test(driver_via_memory driver_via_memory
                    KIND int
                    SOURCE test/Main.cc test/driver_via_memory.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})

    add_amanzi_test(driver_via_memory_restarts driver_via_memory_restarts
                    KIND int
                    SOURCE test/Main.cc test/driver_via_memory_restarts.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})
endif()
