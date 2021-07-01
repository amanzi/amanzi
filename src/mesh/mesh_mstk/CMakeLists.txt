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
project(MESH_MSTK)

# External (TPL) include directories

# Reverse the items in MSTK_INCLUDE_DIRS because when we ask the
# include_directories command to prepend a set of directories using
# the BEFORE keyword, it prepends each item in the MSTK_INCLUDE_DIRS
# list, thereby effectively reversing the order of the
# directories. For example, if we had MYDIRS=a;b;c and we specified
# MYDIRS in the include_directories command, then the
# INCLUDE_DIRECTORIES property for will end up with
# c;b;a;whatever_existed_before

list(REVERSE MESH_MSTK_INCLUDE_DIRS)
include_directories(BEFORE ${MESH_MSTK_INCLUDE_DIRS})

include_directories(${Epetra_INCLUDE_DIRS})
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
# Library: mesh_mstk
#
set(mesh_mstk_files Mesh_MSTK.cc Mesh_MSTK_deform.cc MeshColumn.cc MeshSurfaceCell.cc) 
file(GLOB mesh_mstk_inc_files "*.hh")
add_amanzi_library(mesh_mstk
                   SOURCE ${mesh_mstk_files} HEADERS ${mesh_mstk_inc_files}
                   LINK_LIBS geometry mesh atk ${MSTK_LIBRARIES} ${Epetra_LIBRARIES}
		             ${Teuchos_LIBRARIES} ${Zoltan_LIBRARIES} ${HDF5_LIBRARIES})

if (BUILD_TESTS)
    
    include_directories(${UnitTest_INCLUDE_DIRS})

    # Need to copy files for the tests 
    if (NOT (${MESH_MSTK_SOURCE_DIR} EQUAL ${MESH_MSTK_BINARY_DIR}))
        file(GLOB DataFiles "${MESH_MSTK_SOURCE_DIR}/test/*.xml"
                            "${MESH_MSTK_SOURCE_DIR}/test/*.txt"
                            "${MESH_MSTK_SOURCE_DIR}/test/*.exo"
                            "${MESH_MSTK_SOURCE_DIR}/test/*.par.*")
        file(COPY ${DataFiles} DESTINATION ${MESH_MSTK_BINARY_DIR}/test/)
    endif()

    #Test: eliminate degeneracies for serial mesh
    add_amanzi_test(mesh_mstk_pinchout_test mesh_mstk_pinchout_test
                    KIND int
                    SOURCE test/Main.cc test/mesh_mstk_pinchouts.cc
                    LINK_LIBS mesh_mstk mesh_factory mesh_audit ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})

    #Test: extracted surface from cell from volumetric face label
    add_amanzi_test(mesh_mstk_surface_face_inheritance mesh_mstk_surface_face_inheritance
                    KIND unit
                    SOURCE test/Main.cc test/mesh_mstk_label_inheritance.cc
                    LINK_LIBS mesh_mstk geometry ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})

    #Test: eliminate degeneracies in parallel mesh created from inline partitioning
    add_amanzi_test(mesh_mstk_pinchout_par mesh_mstk_pinchout_par
                    KIND int
		    NPROCS 4
                    SOURCE test/Main.cc test/mesh_mstk_pinchouts.cc
                    LINK_LIBS mesh_mstk mesh_factory mesh_audit ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})

    #Test: eliminate degeneracies in pre-partitioned parallel mesh
    add_amanzi_test(mesh_mstk_pinchout_prepart mesh_mstk_pinchout_prepart
                    KIND int
		    NPROCS 4
                    SOURCE test/Main.cc test/mesh_mstk_degen_prepart.cc
                    LINK_LIBS mesh_mstk mesh_factory mesh_audit ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})

    # Test: mesh_mstk
    add_amanzi_test(mesh_mstk_basis mesh_mstk_basis
                    KIND int
                    SOURCE
                         test/Main.cc 
                         test/mesh_mstk_hex_1x1x1.cc
                         test/mesh_mstk_hex_2x2x1.cc
                         test/mesh_mstk_hex_3x3x3.cc
                         test/mesh_mstk_hex_3x3x3_sets.cc
                         test/mesh_mstk_hex_gen_3x3x3.cc
                         test/mesh_mstk_quad_gen_3x3.cc
                         test/mesh_mstk_extract_surface.cc
 			 test/mesh_mstk_extract_column.cc
                         test/mesh_mstk_deform_vols.cc
		         test/mesh_mstk_extnode_map.cc
                    LINK_LIBS mesh_mstk mesh_audit ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES}) 


    # Test: mesh_mstk_parallel
    add_amanzi_test(mesh_mstk_parallel mesh_mstk_parallel
                    KIND unit
                    NPROCS 4
                    SOURCE
                         test/Main.cc
                         test/mesh_mstk_hex_3x3x3_4P.cc
                         test/mesh_mstk_hex_3x3x3_sets_4P.cc
                         test/mesh_mstk_hex_3x3x3_par_read_4P.cc
                         test/mesh_mstk_quad_gen_3x3_4P.cc
                         test/mesh_mstk_hex_gen_3x3x3_4P.cc
                         test/mesh_mstk_edges_4P.cc
                         test/mesh_mstk_extract_surface_4P.cc
		         test/mesh_mstk_extface_map_4P.cc
		         test/mesh_mstk_extnode_map.cc
                         test/mesh_mstk_deform_vols.cc
                    LINK_LIBS mesh_mstk mesh_audit ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})

    # Test: column_mesh
    add_amanzi_test(mesh_mstk_column mesh_mstk_column
                    KIND unit
                    SOURCE
                          test/Main.cc 
                          test/mesh_mstk_column.cc 
			  LINK_LIBS mesh_mstk geometry ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})

    # Test: extracted surface from column_mesh
    add_amanzi_test(mesh_mstk_surface_column mesh_mstk_surface_column
                    KIND unit
                    SOURCE
                          test/Main.cc 
                          test/mesh_mstk_surface_column.cc 
                    LINK_LIBS mesh_mstk geometry ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})

    # Test: extract fracture test
    add_amanzi_test(mesh_mstk_fracture mesh_mstk_fracture
                    KIND unit
                    SOURCE
                          test/Main.cc 
                          test/mesh_mstk_fracture.cc 
                    LINK_LIBS mesh_mstk geometry ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})
endif()

