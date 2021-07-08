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
include_directories(${ATK_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})

# Library: simple_mesh
file(GLOB mesh_simple_source_files "*.cc")
file(GLOB mesh_simple_inc_files "*.hh")
add_amanzi_library(mesh_simple
                   SOURCE ${mesh_simple_source_files} HEADERS ${mesh_simple_inc_files}
                   LINK_LIBS mesh geometry error_handling atk 
                             ${Epetra_LIBRARIES} ${Teuchos_LIBRARIES})

if (BUILD_TESTS)

    # Add UnitTest include directories
    include_directories(${UnitTest_INCLUDE_DIRS})

    # Test uses Epetra
    include_directories(${Epetra_INCLUDE_DIRS})

   # Need to copy files for the tests 
    if (NOT (${MESH_SIMPLE_SOURCE_DIR} EQUAL ${MESH_SIMPLE_BINARY_DIR}))
      file(GLOB DataFiles "${MESH_SIMPLE_SOURCE_DIR}/test/*.xml"
                          "${MESH_SIMPLE_SOURCE_DIR}/test/*.txt")
      file(COPY ${DataFiles} DESTINATION ${MESH_SIMPLE_BINARY_DIR}/test/)
    endif()

    # Test: mesh_simple
    add_amanzi_test(mesh_simple_all mesh_simple_all
                    KIND unit
		    SOURCE
                          test/Main.cc 
                          test/mesh_simple_maps.cc 
                          test/mesh_simple_cell_numbering.cc 
                          test/mesh_simple_face_adj_cells.cc 
                          test/mesh_simple_sets.cc
                    LINK_LIBS mesh_simple ${UnitTest_LIBRARIES})

endif()

