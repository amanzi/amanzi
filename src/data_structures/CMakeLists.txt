# -*- mode: cmake -*-

#
#  Amanzi
#    Data Structures
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)
include(LibraryManager)

#
# Define a project name
# After this command the following varaibles are defined
#   DATA_STRUCTURES_SOURCE_DIR
#   DATA_STRUCTURES_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(DATA_STRUCTURES)

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${MESH_MSTK_SOURCE_DIR})
include_directories(${MESH_SIMPLE_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${WHETSTONE_SOURCE_DIR})
include_directories(${DATA_STRUCTURES_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${Thyra_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${MSTK_INCLUDE_DIRS})

# Need this define. Errors from MSTK include files 
# about MPI_COMM_WORLD. --lpritch
add_definitions("-DMSTK_HAVE_MPI")


# Library: data_structures
#
#file(GLOB data_structures_src_files "*.cc")
set(data_structures_src_files BlockVector.cc
                              CompositeVectorSpace.cc
                              CompositeVector.cc
                              NoxVector.cc
                              TreeVectorSpace.cc
                              TreeVector.cc
                              ParallelCommunication.cc
                              GraphFE.cc
                              MatrixFE.cc
                              SuperMapLumped.cc
                              SuperMap.cc)

file(GLOB data_structures_inc_files "*.hh")

add_amanzi_library(data_structures 
                   SOURCE ${data_structures_src_files}
                   HEADERS ${data_structures_inc_files}
		   LINK_LIBS mesh error_handling 
                             ${Teuchos_Libraries} ${Epetra_Libraries}
                             ${Thyra_Libraries} ${NOX_LIBRARIES})

#
# Install Target
#

if (BUILD_TESTS)
    # Add UnitTest includes
    include_directories(${UnitTest_INCLUDE_DIRS})

    # Add Data_structures include directory.
    include_directories(${DATA_STRUCTURES_SOURCE_DIR})

    include_directories(${MESH_SOURCE_DIR})
    include_directories(${MESH_SIMPLE_SOURCE_DIR})
    include_directories(${MESH_FACTORY_SOURCE_DIR})

    # Copy test subdirectory for out of source builds
    if (NOT (DATA_STRUCTURES_SOURCE_DIR STREQUAL DATA_STRUCTURES_BINARY_DIR))
        file(GLOB DataFiles "${DATA_STRUCTURES_SOURCE_DIR}/test/*.xml"
                            "${DATA_STRUCTURES_SOURCE_DIR}/test/*.exo")
        file(COPY ${DataFiles} DESTINATION ${DATA_STRUCTURES_BINARY_DIR}/test/)
    endif()

    set(amanzi_libs data_structures geometry mesh_mstk mesh_factory mesh)

    add_amanzi_test(data_structures_vectors data_structures_vectors
                    KIND unit
                    NPROCS 2
		    SOURCE test/Main.cc 
			   test/data_structures_composite_vector.cc
			   test/data_structures_tree_vector.cc
                           test/data_structures_vandelay.cc
                           # test/test_thyra_wrappers.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${NOX_LIBRARIES})

    # Test: tests for FE graph and matrices
    add_amanzi_test(data_structures_graph_fe data_structures_graph_fe
                    KIND unit
                    SOURCE test/Main.cc test/data_structures_graph_fe.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${NOX_LIBRARIES})

    add_amanzi_test(data_structures_matrix_fe data_structures_matrix_fe
                    KIND unit
                    SOURCE test/Main.cc test/data_structures_matrix_fe.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${NOX_LIBRARIES})

    add_amanzi_test(data_structures_matrix_fe_np2 data_structures_matrix_fe NPROCS 2 KIND unit)
    add_amanzi_test(data_structures_matrix_fe_np4 data_structures_matrix_fe NPROCS 4 KIND unit)

    add_amanzi_test(data_structures_supermap data_structures_supermap 
                    KIND unit
                    SOURCE test/Main.cc test/data_structures_supermap.cc test/data_structures_supermap_block.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${NOX_LIBRARIES})
    add_amanzi_test(data_structures_supermap_np3 data_structures_supermap NPROCS 3 KIND unit)
    add_amanzi_test(data_structures_supermap_np5 data_structures_supermap NPROCS 5 KIND unit)

    add_amanzi_test(data_structures_supermap_block data_structures_supermap_block
                    KIND unit
                    SOURCE test/Main.cc test/data_structures_supermap_block.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${NOX_LIBRARIES})

    # add_amanzi_test(data_structures_matrix_fe_speed data_structures_matrix_fe_speed
    #                 KIND unit
    # 	              SOURCE test/Main.cc test/data_structures_matrix_fe_speed.cc
    # 	              LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    # add_amanzi_test(data_structures_access data_structures_access_times
    #                 KIND int
    #		      SOURCE test/Main.cc test/data_structures_cv_access.cc
    #		      LINK_LIBS data_structures mesh
    #	                        geometry mesh_factory
    #                           simple_mesh
    #	                        ${UnitTest_LIBRARIES})
endif()
