# -*- mode: cmake -*-

#
#  Amanzi
#    MSTK Mesh frameworks
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)

#
# Define a project name
# After this command the following varaibles are defined
#   MSTK_SOURCE_DIR
#   MSTK_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(MSTK)

# External (TPL) include directories

# Reverse the items in MSTK_INCLUDE_DIRS because when we ask the
# include_directories command to prepend a set of directories using
# the BEFORE keyword, it prepends each item in the MSTK_INCLUDE_DIRS
# list, thereby effectively reversing the order of the
# directories. For example, if we had MYDIRS=a;b;c and we specified
# MYDIRS in the include_directories command, then the
# INCLUDE_DIRECTORIES property for will end up with
# c;b;a;whatever_existed_before

list(REVERSE MSTK_INCLUDE_DIRS)
include_directories(BEFORE ${MSTK_INCLUDE_DIRS})

include_directories(${Tpetra_INCLUDE_DIRS})
include_directories(${Teuchos_INCLUDE_DIRS})

# Amanzi include directories
include_directories(${ATK_SOURCE_DIR})
include_directories(${DBC_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})


# Need this define. Errors from MSTK include files
# about MPI_COMM_WORLD. --lpritch
add_definitions("-DMSTK_HAVE_MPI")

if (WITH_MSTK_2_20rc1_OR_NEWER)
  add_definitions("-DMSTK_2_20rc1_OR_NEWER")
endif ()

if (WITH_MSTK_2_21rc1_OR_NEWER)
  add_definitions("-DMSTK_2_21rc1_OR_NEWER")
endif ()

if (WITH_MSTK_3_00_OR_NEWER)
  add_definitions("-DMSTK_3_00_OR_NEWER")
endif ()

#
# Library: mstk_mesh
#
set(mstk_mesh_files Mesh_MSTK.cc Mesh_MSTK_deform.cc MeshColumn.cc MeshSurfaceCell.cc)
file(GLOB mstk_inc_files "*.hh")
add_amanzi_library(mstk_mesh
                   SOURCE ${mstk_mesh_files} HEADERS ${mstk_inc_files}
                   LINK_LIBS geometry mesh atk ${MSTK_LIBRARIES} ${Tpetra_LIBRARIES}
		             ${Teuchos_LIBRARIES} ${Zoltan_LIBRARIES} ${HDF5_LIBRARIES})

if (BUILD_TESTS)

    include_directories(${UnitTest_INCLUDE_DIRS})

    # Need to copy files for the tests
    if (NOT (${MSTK_SOURCE_DIR} EQUAL ${MSTK_BINARY_DIR}))
        file(GLOB DataFiles "${MSTK_SOURCE_DIR}/test/*.xml"
                            "${MSTK_SOURCE_DIR}/test/*.txt"
                            "${MSTK_SOURCE_DIR}/test/*.exo"
                            "${MSTK_SOURCE_DIR}/test/*.par.*")
        file(COPY ${DataFiles} DESTINATION ${MSTK_BINARY_DIR}/test/)
    endif()

    #Test: eliminate degeneracies for serial mesh
    # add_amanzi_test(mstk_pinchout_test mstk_pinchout_test
    #                 KIND int
    #                 SOURCE test/Main.cc 
    #                 test/test_pinchouts.cc
    #                 LINK_LIBS mstk_mesh mesh_factory mesh_audit ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})

    # #Test: eliminate degeneracies in parallel mesh created from inline partitioning
    # add_amanzi_test(mstk_pinchout_test_par mstk_pinchout_test_par
    #                 --kokkos-ndevices=1
    #     	    KIND int
    #     	    NPROCS 4
    #                 SOURCE test/Main.cc test/test_pinchouts.cc
    #                 LINK_LIBS mstk_mesh mesh_factory mesh_audit ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})

    #Test: eliminate degeneracies in pre-partitioned parallel mesh
    # add_amanzi_test(mstk_pinchout_test_prepart mstk_pinchout_test_prepart
    #     	    --kokkos-ndevices=1
    #                 KIND int
    #     	    NPROCS 4
    #                 SOURCE test/Main.cc test/test_degen_prepart.cc
    #                 LINK_LIBS mstk_mesh mesh_factory mesh_audit ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})

    # Test: mstk_mesh
    add_amanzi_test(mstk_mesh_basis mstk_mesh_basis
                    KIND int
                    SOURCE
                         test/Main.cc
                         test/test_hex_1x1x1.cc
                         test/test_hex_2x2x1.cc
                         test/test_hex_3x3x3.cc
                         test/test_hex_3x3x3_sets.cc
                         test/test_hex_gen_3x3x3.cc
                         test/test_quad_gen_3x3.cc
                         test/test_extract_surface.cc
			 test/test_extract_column.cc
#                         test/test_deform_vols.cc
                    LINK_LIBS mstk_mesh mesh_audit ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})


    # Test: mstk_mesh_parallel
    add_amanzi_test(mstk_mesh_parallel mstk_mesh_parallel
		    --kokkos-ndevices=1
                    KIND unit
                    NPROCS 4
                    SOURCE
                         test/Main.cc
                         test/test_hex_3x3x3_4P.cc
                         test/test_hex_3x3x3_sets_4P.cc
                         test/test_hex_3x3x3_par_read_4P.cc
                         test/test_quad_gen_3x3_4P.cc
                         test/test_hex_gen_3x3x3_4P.cc
                         test/test_edges_4P.cc
                         test/test_extract_surface_4P.cc
		         test/test_extface_map_4P.cc
                         test/test_deform_vols.cc
                    LINK_LIBS mstk_mesh mesh_audit ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})

    # Test: column_mesh
    add_amanzi_test(mstk_column_mesh mstk_column_mesh
                    KIND unit
                    SOURCE
                          test/Main.cc
                          test/test_column_mesh.cc
			  LINK_LIBS mstk_mesh geometry ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})
    # Test: extracted surface from column_mesh
    add_amanzi_test(mstk_surface_column_mesh mstk_surface_column_mesh
                    KIND unit
                    SOURCE
                          test/Main.cc
                          test/test_surface_column_mesh.cc
                    LINK_LIBS mstk_mesh geometry ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})

endif()
