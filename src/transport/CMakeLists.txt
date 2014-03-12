# -*- mode: cmake -*-

#
#  Amanzi
#   Transport process kernel
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)

#
# Define a project name
# After this command the following varaibles are defined
#   TRANSPORT_SOURCE_DIR
#   TRANSPORT_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(TRANSPORT)

# Get the includes directory for chemistry (needed for Chemistry_State and BCs).
get_property(CHEM_INCLUDES_DIR GLOBAL PROPERTY CHEM_INCLUDES_DIR)

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${MESH_DATA_SOURCE_DIR})
include_directories(${MESH_FACTORY_SOURCE_DIR})
include_directories(${DATA_STRUCTURES_SOURCE_DIR})
include_directories(${STATE_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${WHETSTONE_SOURCE_DIR})
include_directories(${SOLVERS_SOURCE_DIR})
include_directories(${OUTPUT_SOURCE_DIR})
include_directories(${CHEM_INCLUDES_DIR})
include_directories(${CHEMPK_SOURCE_DIR})
include_directories(${TIME_INTEGRATION_SOURCE_DIR})
include_directories(${FUNCS_SOURCE_DIR})
include_directories(${MFUNCS_SOURCE_DIR})
include_directories(${DBG_SOURCE_DIR})
include_directories(${FLOW_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${ASCEMIO_INCLUDE_DIR})

#
# Library: transport
#
set(transport_inc_files
  Dispersion.hh
  Dispersion_TPFA.hh
  Reconstruction.hh
  TransportBCFactory.hh
  TransportBoundaryFunction.hh
  TransportBoundaryFunction_Tracer.hh
  TransportBoundaryFunction_Alquimia.hh
  TransportDefs.hh
  TransportDomainFunction.hh
  Transport_PK.hh
  TransportSourceFactory.hh)

set(transport_src_files Transport_PK.cc Transport_TI.cc 
                        Transport_Limiters.cc Transport_VandV.cc Transport_IO.cc 
                        Reconstruction.cc 
                        DispersionMatrixFactory.cc Dispersion.cc 
                        Dispersion_TPFA.cc Dispersion_MFD.cc Dispersion_NLFV.cc
                        TransportBCFactory.cc TransportSourceFactory.cc
                        TransportBoundaryFunction_Tracer.cc TransportBoundaryFunction_Alquimia.cc
                        TransportDomainFunction.cc)

add_amanzi_library(transport 
                   SOURCE ${transport_src_files} HEADERS ${transport_inc_files}
                   LINK_LIBS geometry mesh error_handling state data_structures 
                             geochemistry chemistry_pk time_integration ${Epetra_LIBRARIES})

if ( BUILD_TESTS )

    # Add UnitTest includes
    include_directories(${UnitTest_INCLUDE_DIRS})

    if (NOT (${TRANSPORT_SOURCE_DIR} EQUAL ${TRANSPORT_BINARY_DIR}) )
        execute_process(COMMAND ${CMAKE_COMMAND} -E 
          copy_directory ${TRANSPORT_SOURCE_DIR}/test ${TRANSPORT_BINARY_DIR}/test) 
    endif()

    include_directories(${TRANSPORT_SOURCE_DIR})
    include_directories(${Boost_INCLUDE_DIRS})

    set(amanzi_libs transport state mesh_functions data_structures functions whetstone mesh_factory mesh_audit mesh geometry)

    # Test: miscaleneous transport routines
    add_amanzi_test(transport_mics transport_mics
                    KIND unit
                    SOURCE test/Main.cc test/transport_mics.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    add_amanzi_test(transport_advance_simple transport_advance_simple
                    KIND unit
                    SOURCE test/Main.cc test/transport_advance_simple.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    add_amanzi_test(transport_limiters transport_limiters
                    KIND unit
                    SOURCE test/Main.cc test/transport_limiters.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    add_amanzi_test(transport_convergence transport_convergence
                    KIND unit
                    SOURCE test/Main.cc test/transport_convergence.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    add_amanzi_test(transport_dispersion transport_dispersion
                    KIND unit
                    SOURCE test/Main.cc test/transport_dispersion.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    if ( ENABLE_MSTK_Mesh )
        set(amanzi_libs geometry mesh mesh_audit mesh_factory state transport whetstone output mesh_functions functions data_structures)

        add_amanzi_test(transport_advance_mstk transport_advance_mstk
	                KIND unit
                        SOURCE test/Main.cc test/transport_advance_mstk.cc
                        LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

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

        add_amanzi_test(transport_matrix transport_matrix
	                KIND unit
                        SOURCE test/Main.cc test/transport_matrix.cc
                        LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})
    endif()    

    if ( ENABLE_STK_Mesh )
        set(amanzi_libs geometry mesh mesh_audit mesh_factory state transport whetstone mesh_functions data_structures)

        add_amanzi_test(transport_advance_stk transport_advance_stk
	                KIND unit
                        SOURCE test/Main.cc test/transport_advance_stk.cc
                        LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})
        
        add_amanzi_test(transport_parallel_stk transport_parallel_stk
			KIND unit
	                NPROCS 4
                        SOURCE test/Main.cc test/transport_parallel_stk.cc
                        LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})
          
    endif()    
endif()


