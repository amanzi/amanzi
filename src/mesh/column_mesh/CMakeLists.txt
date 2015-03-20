# -*- mode: cmake -*-

#
#  Amanzi
#    Column Mesh Class
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)


#
# Define a project name
# After this command the following varaibles are defined
#   COLUMN_MESH_SOURCE_DIR
#   COLUMN_MESH_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(COLUMN_MESH)

# Amanzi include directories
include_directories(${MESH_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})

# Library: column_mesh
file(GLOB column_mesh_source_files "*.cc")
file(GLOB column_mesh_inc_files "*.hh")
add_amanzi_library(column_mesh
                   SOURCE ${column_mesh_source_files} HEADERS ${column_mesh_inc_files}
                   LINK_LIBS mesh geometry error_handling ${Epetra_LIBRARIES} ${Teuchos_LIBRARIES})

if (BUILD_TESTS)

    # Add UnitTest include directories
    include_directories(${UnitTest_INCLUDE_DIRS})

    # Need mesh factory
    include_directories(${MESH_FACTORY_SOURCE_DIR})

   # Need to copy files for the tests 
    if (NOT (${COLUMN_MESH_SOURCE_DIR} EQUAL ${COLUMN_MESH_BINARY_DIR}))
        execute_process(COMMAND ${CMAKE_COMMAND} -E 
          copy_directory ${COLUMN_MESH_SOURCE_DIR}/test ${COLUMN_MESH_BINARY_DIR}/test) 
    endif()

    # Test: column_mesh
    add_amanzi_test(column_mesh test_column_mesh
                    KIND unit
		    SOURCE
                          test/Main.cc 
                          test/test_column_mesh.cc 
                    LINK_LIBS column_mesh mesh_factory geometry ${UnitTest_LIBRARIES})

endif()

