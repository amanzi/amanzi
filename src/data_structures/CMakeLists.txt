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
include_directories(${MSTK_SOURCE_DIR})
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
set(data_structures_src_files BlockSpace.cc
  CompositeSpace.cc
  CompositeVectorSpace.cc
  CompositeVector.cc
  #NoxVector.cc
  TreeVectorSpace.cc
  TreeVector.cc
  GraphFE.cc
  MatrixFE.cc
  SuperMapLumped.cc
  SuperMap.cc
  )

set(data_structures_inc_files DataStructuresHelpers.hh
  VectorHarness.hh
  CSR.hh
  BlockSpace.hh
  BlockVector_decl.hh
  BlockVector_impl.hh
  BlockVector.hh
  CompositeSpace.hh
  CompositeVector_decl.hh
  CompositeVector_impl.hh
  CompositeVector.hh 
  CompositeVectorSpace.hh
  CompositeMatrix.hh 
  TreeVectorSpace.hh
  TreeVector_Utils.hh
  TreeVector.hh
  TreeMatrix.hh
  TensorVector.hh
  GraphFE.hh
  MatrixFE.hh
  #NoxVector.hh
  SuperMap.hh
  SuperMapLumped.hh
  Patch.hh
  )

set(data_structures_link_libs mesh
  error_handling 
  ${Teuchos_Libraries}
  ${Tpetra_Libraries}
  #${Thyra_Libraries}
  #${NOX_LIBRARIES}
  )


#
# Install Target
#
add_amanzi_library(data_structures 
                   SOURCE ${data_structures_src_files}
                   HEADERS ${data_structures_inc_files}
		   LINK_LIBS ${data_structures_link_libs}
                   )


if (BUILD_TESTS)
    # Add UnitTest includes
    include_directories(${UnitTest_INCLUDE_DIRS})

    # Add Data_structures include directory.
    include_directories(${DATA_STRUCTURES_SOURCE_DIR})

    include_directories(${MESH_SOURCE_DIR})
    include_directories(${MESH_SIMPLE_SOURCE_DIR})
    include_directories(${MESH_FACTORY_SOURCE_DIR})
    include_directories(${MSTK_SOURCE_DIR})

    # Copy test subdirectory for out of source builds
    if (NOT (DATA_STRUCTURES_SOURCE_DIR STREQUAL DATA_STRUCTURES_BINARY_DIR))
        file(GLOB DataFiles "${DATA_STRUCTURES_SOURCE_DIR}/test/*.xml"
                            "${DATA_STRUCTURES_SOURCE_DIR}/test/*.exo")
        file(COPY ${DataFiles} DESTINATION ${DATA_STRUCTURES_BINARY_DIR}/test/)
    endif()

    set(data_structures_test_link_libs data_structures geometry mstk_mesh mesh_factory mesh)

    add_amanzi_test(data_structures_vectors data_structures_vectors
    		    --kokkos-ndevices=1
                    KIND unit
                    NPROCS 2
        	    SOURCE test/Main.cc
                           test/data_structures_communication.cc
        		   test/data_structures_block_vector.cc
        		   test/data_structures_composite_vector.cc
        		   test/data_structures_int_composite_vector.cc
        		   test/data_structures_tree_vector.cc
                           test/data_structures_vandelay.cc
                    LINK_LIBS ${data_structures_test_link_libs} ${UnitTest_LIBRARIES})

    # Test: tests for FE graph and matrices
    add_amanzi_test(data_structures_graph_fe data_structures_graph_fe
                    KIND unit
                    SOURCE test/Main.cc test/data_structures_graph_fe.cc
                    LINK_LIBS ${data_structures_test_link_libs} ${UnitTest_LIBRARIES})

    add_amanzi_test(data_structures_matrix_fe data_structures_matrix_fe
                    KIND unit
                    SOURCE test/Main.cc test/data_structures_matrix_fe.cc
                    LINK_LIBS ${data_structures_test_link_libs} ${UnitTest_LIBRARIES})

    add_amanzi_test(data_structures_matrix_fe_np2 data_structures_matrix_fe NPROCS 2 KIND unit)
    add_amanzi_test(data_structures_matrix_fe_np4 data_structures_matrix_fe NPROCS 4 KIND unit)

    add_amanzi_test(data_structures_supermap data_structures_supermap
                    KIND unit
                    SOURCE test/Main.cc test/data_structures_supermap.cc
                    LINK_LIBS ${data_structures_test_link_libs} ${UnitTest_LIBRARIES})
    add_amanzi_test(data_structures_supermap_np3 data_structures_supermap --kokkos-ndevices=1 NPROCS 3 KIND unit)
    add_amanzi_test(data_structures_supermap_np5 data_structures_supermap --kokkos-ndevices=1 NPROCS 5 KIND unit)

    # add_amanzi_test(data_structures_block_supermap data_structures_block_supermap
    #                 KIND unit
    #                 SOURCE test/Main.cc test/data_structures_block_supermap.cc
    #                 LINK_LIBS ${data_structures_test_link_libs} ${UnitTest_LIBRARIES})
    # add_amanzi_test(test_block_supermap_np3 test_block_supermap NPROCS 3 KIND unit)
    # add_amanzi_test(test_block_supermap_np5 test_block_supermap NPROCS 5 KIND unit)

    # add_amanzi_test(matrix_fe_test_speed matrix_fe_test_speed
    #                 KIND unit
    # 	              SOURCE test/Main.cc test/data_structures_matrix_fe_speed.cc
    # 	              LINK_LIBS ${data_structures_test_link_libs} ${UnitTest_LIBRARIES})

    add_amanzi_test(data_structures_common_patterns, data_structures_common_patterns
                    KIND unit
                    SOURCE test/Main.cc test/data_structures_common_patterns.cc
                    LINK_LIBS data_structures mesh geometry mesh_factory simple_mesh ${UnitTest_LIBRARIES})

    add_amanzi_test(data_structures_patch, data_structures_patch
                    KIND unit
                    SOURCE test/Main.cc test/data_structures_patch.cc
                    LINK_LIBS data_structures mesh geometry mesh_factory simple_mesh ${UnitTest_LIBRARIES})
endif()
