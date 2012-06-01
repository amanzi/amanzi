# -*- mode: cmake -*-

#
#  Amanzi
#    Mesh Base Class
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)


#
# Define a project name
# After this command the following varaibles are defined
#   MESH_SIMPLE_SOURCE_DIR
#   MESH_SIMPLE_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(MESH_SIMPLE)

# Amanzi include directories
include_directories(${MESH_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})

# Library: simple_mesh
file(GLOB simple_mesh_source_files "*.cc")
add_library(simple_mesh ${simple_mesh_source_files})
add_install_library(simple_mesh)
target_link_libraries(simple_mesh mesh geometry error_handling ${Epetra_LIBRARIES} ${Teuchos_LIBRARIES})

#
# Install Targets
#
file(GLOB simple_inc_files "*.hh")
add_install_include_file(${simple_inc_files})

if (BUILD_TESTS)

    # Add UnitTest include directories
    include_directories(${UnitTest_INCLUDE_DIRS})

    # Test uses Epetra
    include_directories(${Epetra_INCLUDE_DIRS})

   # Need to copy files for the tests 
    if (NOT (${MESH_SIMPLE_SOURCE_DIR} EQUAL ${MESH_SIMPLE_BINARY_DIR}))
        execute_process(COMMAND ${CMAKE_COMMAND} -E 
          copy_directory ${MESH_SIMPLE_SOURCE_DIR}/test ${MESH_SIMPLE_BINARY_DIR}/test) 
    endif()

    # Test: simple_mesh
    add_amanzi_test(simple_mesh test_simple_mesh
                    KIND unit
		    SOURCE
                          test/Main.cc 
                          test/test_maps.cpp 
                          test/test_cell_numbering.cc 
                          test/test_face_adj_cells.cc 
                          test/test_node_cell_faces.cc
                          test/test_sets.cc
                    LINK_LIBS simple_mesh ${UnitTest_LIBRARIES})

endif()

