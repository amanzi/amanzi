# -*- mode: cmake -*-

#
#  Amanzi
#   Transport process kernel
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
##include(PrintVariable)
##include(TestManager)

#
# Define a project name
# After this command the following varaibles are defined
#   TRANSPORT_SOURCE_DIR
#   TRANSPORT_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
##project(TRANSPORT)

# Get the includes directory for chemistry (needed for BCs).
get_property(CHEM_INCLUDES_DIR GLOBAL PROPERTY CHEM_INCLUDES_DIR)

# External (TPL) include directories
##include_directories(${Epetra_INCLUDE_DIRS})
##include_directories(${ASCEMIO_INCLUDE_DIR})
# Amanzi include directories
##include_directories(${ATK_SOURCE_DIR})
##include_directories(${CHEM_INCLUDES_DIR})
##include_directories(${CHEMPK_SOURCE_DIR})
##include_directories(${DATA_STRUCTURES_SOURCE_DIR})
##include_directories(${FUNCS_SOURCE_DIR})
##include_directories(${GEOMETRY_SOURCE_DIR})
##include_directories(${MESH_SOURCE_DIR})
##include_directories(${MESH_FACTORY_SOURCE_DIR})
##include_directories(${MFUNCS_SOURCE_DIR})
##include_directories(${OUTPUT_SOURCE_DIR})
##include_directories(${OPERATORS_SOURCE_DIR})
##include_directories(${PKS_SOURCE_DIR})
##include_directories(${SOLVERS_SOURCE_DIR})
##include_directories(${STATE_SOURCE_DIR})
##include_directories(${TIME_INTEGRATION_SOURCE_DIR})
##include_directories(${WHETSTONE_SOURCE_DIR})

include_directories(${Amanzi_TPL_MSTK_INCLUDE_DIRS})

#
# Transport registrations
#



#
# Library: transport
#
set(transport_inc_files
  DiffusionPhase.hh
  MaterialProperties.hh
  MDM.hh
  MDMFactory.hh
  MDMPartition.hh
  MDM_Isotropic.hh
  MDM_Bear.hh
  MultiscaleTransportPorosity.hh
  MultiscaleTransportPorosityFactory.hh
  MultiscaleTransportPorosityPartition.hh
  MultiscaleTransportPorosity_DPM.hh
#  TransportBCFactory.hh
  TransportDomainFunction.hh 
  TransportBoundaryFunction_Alquimia.hh
  TransportSourceFunction_Alquimia.hh
  TransportDefs.hh
  Transport_PK_ATS.hh)

set(transport_src_files Transport_PK.cc Transport_TI.cc 
                        Transport_VandV.cc Transport_Initialize.cc 
                        Transport_Dispersion.cc Transport_HenryLaw.cc
                        MDM_Isotropic.cc MDM_Bear.cc MDM_BurnettFrind.cc MDM_LichtnerKelkarRobinson.cc
                        MDMPartition.cc MDMFactory.cc
                        MultiscaleTransportPorosityFactory.cc MultiscaleTransportPorosity_DPM.cc
                        MultiscaleTransportPorosityPartition.cc
                        TransportBoundaryFunction_Alquimia.cc
                        TransportSourceFunction_Alquimia.cc)

#                        TransportBoundaryFunction_Alquimia.cc)
#                        TransportBCFactory.cc #TransportSourceFactory.cc
#                        TransportBoundaryFunction_Tracer.cc TransportBoundaryFunction_Alquimia.cc
#                        TransportBoundaryFunction_Coupler.cc)

if (ENABLE_ALQUIMIA)
  set(geochem_lib geochemistry)
endif()

add_library(pk_transport 
            ${transport_src_files}
)

install(TARGETS pk_transport DESTINATION lib)

#                   LINK_LIBS geometry mesh error_handling state data_structures 
#                             ${geochem_lib} chemistry_pk time_integration ${Epetra_LIBRARIES})

if ( 0 )
#if ( BUILD_TESTS )
    # Add UnitTest includes
    include_directories(${UnitTest_INCLUDE_DIRS})

    # Copy test directory files if any out of source build
    if (NOT (${TRANSPORT_SOURCE_DIR} EQUAL ${TRANSPORT_BINARY_DIR}) )
        # execute_process(COMMAND ${CMAKE_COMMAND} -E 
        #   copy_directory ${TRANSPORT_SOURCE_DIR}/test ${TRANSPORT_BINARY_DIR}/test) 
        file(GLOB DataFiles "${TRANSPORT_SOURCE_DIR}/test/*.xml"
                            "${TRANSPORT_SOURCE_DIR}/test/*.exo"
                            "${TRANSPORT_SOURCE_DIR}/test/*.par.*")
        file(COPY ${DataFiles} DESTINATION ${TRANSPORT_BINARY_DIR}/test/)
    endif()

    include_directories(${TRANSPORT_SOURCE_DIR})
    include_directories(${Boost_INCLUDE_DIRS})

    set(amanzi_libs transport pks state operators mesh_functions data_structures functions 
                    whetstone mesh_factory mesh_audit mesh geometry)

    # Test: miscaleneous transport routines
    add_amanzi_test(transport_mics transport_mics
                    KIND unit
                    SOURCE test/Main.cc test/transport_mics.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    add_amanzi_test(transport_advance transport_advance
                    KIND unit
                    SOURCE test/Main.cc test/transport_advance.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    #add_amanzi_test(transport_limiters transport_limiters
    #                KIND unit
    #                SOURCE test/Main.cc test/transport_limiters.cc
    #                LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    add_amanzi_test(transport_convergence transport_convergence
                    KIND unit
                    SOURCE test/Main.cc test/transport_convergence.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    add_amanzi_test(transport_dispersion transport_dispersion
                    KIND unit
                    SOURCE test/Main.cc test/transport_dispersion.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    add_amanzi_test(transport_2D transport_2D
                    KIND unit
                    SOURCE test/Main.cc test/transport_2D.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    add_amanzi_test(transport_2D_long transport_2D_long 
                    KIND unit
                    SOURCE test/Main.cc test/transport_2D_long.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    add_amanzi_test(transport_2D_source transport_2D_source
                    KIND unit
                    SOURCE test/Main.cc test/transport_2D_source.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    add_amanzi_test(transport_3D_long transport_3D_long
                    KIND unit
                    SOURCE test/Main.cc test/transport_3D_long.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    add_amanzi_test(transport_subcycling transport_subcycling
                    KIND unit
                    SOURCE test/Main.cc test/transport_subcycling.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    if ( ENABLE_MSTK_Mesh )
	add_amanzi_test(transport_parallel_mstk transport_parallel_mstk
			       KIND int
	               NPROCS 4 
	               SOURCE test/Main.cc test/transport_parallel_mstk.cc
	               LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

        add_amanzi_test(transport_parallel_read_mstk transport_parallel_read_mstk
		       KIND int
	               NPROCS 4 
                       SOURCE test/Main.cc test/transport_parallel_read_mstk.cc
                       LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})
    endif()    
endif()

