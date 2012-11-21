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

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})


# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${MSTK_INCLUDE_DIRS})

# Need this define. Errors from MSTK include files 
# about MPI_COMM_WORLD. --lpritch
add_definitions("-DMSTK_HAVE_MPI")

#
# Library: mstk_mesh
#
set(mstk_mesh_files Mesh_MSTK.cc) 
file(GLOB mstk_inc_files "*.hh")
add_amanzi_library(mstk_mesh
                   SOURCE ${mstk_mesh_files} HEADERS ${mstk_inc_files}
                   LINK_LIBS geometry mesh ${MSTK_LIBRARIES} ${Epetra_LIBRARIES} ${Teuchos_LIBRARIES} ${Zoltan_LIBRARIES})

if (BUILD_TESTS)
    
    # Add UnitTest ExodusII include directories
    include_directories(${UnitTest_INCLUDE_DIRS})

    # Need to copy files for the tests 
    if (NOT (${MSTK_SOURCE_DIR} EQUAL ${MSTK_BINARY_DIR}))
        execute_process(COMMAND ${CMAKE_COMMAND} -E 
          copy_directory ${MSTK_SOURCE_DIR}/test ${MSTK_BINARY_DIR}/test) 
    endif()


    # Test: mstk_mesh
    add_amanzi_test(mstk_mesh test_mstk_mesh
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
                    LINK_LIBS mstk_mesh mesh_audit ${UnitTest_LIBRARIES}) 

    # Test: mstk_mesh_parallel
    add_amanzi_test(mstk_mesh_parallel test_parallel_mstk_mesh
                    KIND unit
                    NPROCS 4
                    SOURCE
                         test/Main.cc
                         test/test_hex_3x3x3_4P.cc
                         test/test_hex_3x3x3_sets_4P.cc
			 test/test_hex_3x3x3_par_read_4P.cc
                         test/test_quad_gen_3x3_4P.cc
                         test/test_hex_gen_3x3x3_4P.cc
			 test/test_extract_surface_4P.cc
                    LINK_LIBS mstk_mesh mesh_audit ${UnitTest_LIBRARIES})

endif()

