# -*- mode: cmake -*-

#
#  Amanzi
#    Interoperability
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)
include(LibraryManager)

#
# Define a project name
# After this command the following varaibles are defined
#   INTEROPERABILITY_SOURCE_DIR
#   INTEROPERABILITY_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(INTEROPERABILITY)

# Amanzi include directories
include_directories(${ATK_SOURCE_DIR})
include_directories(${CHEM_SOURCE_DIR})
include_directories(${CHEMPK_SOURCE_DIR})
include_directories(${DATA_STRUCTURES_SOURCE_DIR})
include_directories(${DBC_SOURCE_DIR})
include_directories(${DBG_SOURCE_DIR})
include_directories(${EOS_SOURCE_DIR})
include_directories(${FUNCS_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${GEOCHEM_SOURCE_DIR})
include_directories(${MFUNCS_SOURCE_DIR})
include_directories(${MESH_DATA_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${MPC_SOURCE_DIR})
include_directories(${MPC_TREE_SOURCE_DIR})
include_directories(${MULTIPHASE_SOURCE_DIR})
include_directories(${OPERATORS_SOURCE_DIR})
include_directories(${OUTPUT_SOURCE_DIR})
include_directories(${PKS_SOURCE_DIR})
include_directories(${SOLVERS_SOURCE_DIR})
include_directories(${STATE_SOURCE_DIR})
include_directories(${STATE_SOURCE_DIR}/data)
include_directories(${STATE_SOURCE_DIR}/evaluators)
include_directories(${TIME_INTEGRATION_SOURCE_DIR})
include_directories(${TRANSPORT_SOURCE_DIR})
include_directories(${WHETSTONE_SOURCE_DIR})

get_property(CHEM_INCLUDES_DIR GLOBAL PROPERTY CHEM_INCLUDES_DIR)
include_directories(${CHEM_INCLUDES_DIR})

# External (TPL) include directories
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${HDF5_C_INCLUDE_DIR})
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${ASCEMIO_INCLUDE_DIR})

# Interoperability directories
include_directories(${ATS_SOURCE_DIR}/src/constitutive_relations)
include_directories(${ATS_SOURCE_DIR}/src/constitutive_relations/column_integrators)
include_directories(${ATS_SOURCE_DIR}/src/constitutive_relations/generic_evaluators)
include_directories(${ATS_SOURCE_DIR}/src/constitutive_relations/eos)
include_directories(${ATS_SOURCE_DIR}/src/constitutive_relations/surface_subsurface_fluxes)
include_directories(${ATS_SOURCE_DIR}/src/operators/upwinding)
include_directories(${ATS_SOURCE_DIR}/src/pks)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/elevation)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/porosity)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/overland_conductivity)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/sources)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/water_content)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/wrm)

if (BUILD_TESTS)
    include_directories(${MESH_FACTORY_SOURCE_DIR})

    # Add UnitTest includes
    include_directories(${UnitTest_INCLUDE_DIRS})
    # Add MPC include directory. Need to remove this required path. -- lpritch
    include_directories(${INTEROPERABILITY_SOURCE_DIR})

    # Copy test subdirectory for out of source builds
    if (NOT (INTEROPERABILITY_SOURCE_DIR STREQUAL INTEROPERABILITY_BINARY_DIR))
      file(GLOB DataFiles "${INTEROPERABILITY_SOURCE_DIR}/test/*.xml"
                          "${INTEROPERABILITY_SOURCE_DIR}/test/*.in"
                          "${INTEROPERABILITY_SOURCE_DIR}/test/*.dat"
                          "${INTEROPERABILITY_SOURCE_DIR}/test/*.exo"
                          "${INTEROPERABILITY_SOURCE_DIR}/test/*.h5")
      file(COPY ${DataFiles} DESTINATION ${INTEROPERABILITY_BINARY_DIR}/test/)
    endif()

    set(amanzi_libs mpc pks mpc_tree transport  
                    pks eos state operators mesh_extracted mesh_factory mesh_audit 
                    geometry whetstone mesh_functions functions)

    set(ats_libs ats_flow ats_eos ats_column_integrator ats_generic_evals ats_surf_subsurf)

    # Test: mpc driver for flow and transport
    add_amanzi_test(interoperability_flow_transport interoperability_flow_transport
                    KIND      int
                    SOURCE    test/Main.cc test/interoperability_flow_transport.cc
                    LINK_LIBS ${amanzi_libs} ${ats_libs}
                              ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES} ${HYPRE_LIBRARY})

endif()
