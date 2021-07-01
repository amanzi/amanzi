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

    add_amanzi_test(data_structures test_data_structures 
                    KIND unit
                    NPROCS 2
		    SOURCE test/Main.cc 
			   test/test_composite_vector.cc
			   test/test_tree_vector.cc
                           test/test_vandelay.cc
                           # test/test_thyra_wrappers.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${NOX_LIBRARIES})

    # Test: tests for FE graph and matrices
    add_amanzi_test(graph_fe_test graph_fe_test
                    KIND unit
                    SOURCE test/Main.cc test/graph_fe_test.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${NOX_LIBRARIES})

    add_amanzi_test(matrix_fe_test matrix_fe_test
                    KIND unit
                    SOURCE test/Main.cc test/matrix_fe_test.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${NOX_LIBRARIES})

    add_amanzi_test(matrix_fe_test_np2 matrix_fe_test NPROCS 2 KIND unit)
    add_amanzi_test(matrix_fe_test_np4 matrix_fe_test NPROCS 4 KIND unit)

    add_amanzi_test(supermap_test supermap_test 
                    KIND unit
                    SOURCE test/Main.cc test/supermap_test.cc test/supermap_block_test.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${NOX_LIBRARIES})
    add_amanzi_test(supermap_test_np3 supermap_test NPROCS 3 KIND unit)
    add_amanzi_test(supermap_test_np5 supermap_test NPROCS 5 KIND unit)

    add_amanzi_test(supermap_block_test supermap_block_test
                    KIND unit
                    SOURCE test/Main.cc test/supermap_block_test.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${NOX_LIBRARIES})

    # add_amanzi_test(matrix_fe_test_speed matrix_fe_test_speed
    #                 KIND unit
    # 	              SOURCE test/Main.cc test/matrix_fe_test_speed.cc
    # 	              LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    # add_amanzi_test(test_access test_access_times
    #                 KIND int
    #		      SOURCE test/Main.cc test/test_cv_access.cc
    #		      LINK_LIBS data_structures mesh
    #	                        geometry mesh_factory
    #                           simple_mesh
    #	                        ${UnitTest_LIBRARIES})
endif()
