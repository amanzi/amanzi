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
target_link_libraries(simple_mesh mesh geometry error_handling ${Epetra_LIBRARIES} ${Teuchos_LIBRARIES})
link_list_add(simple_mesh)
#
# Install Targets
#
file(GLOB simple_inc_files "*.hh")
add_install_include_file(${simple_inc_files})
add_install_library(simple_mesh)

if (BUILD_TESTS)

    # Add UnitTest include directories
    include_directories(${UnitTest_INCLUDE_DIRS})

    # Test uses Epetra
    include_directories(${Epetra_INCLUDE_DIRS})

    # Test: simple_mesh
    # Test does not build. Errors about State. This is a class
    # defined in the mpc and that is not part of the build.
    # Will comment out for now and create a ticket #329
    # to remind me to change this when the mpc is back in the build system
    add_executable(test_simple_mesh test/Main.cc test/test_maps.cpp test/test_cell_numbering.cc test/test_face_adj_cells.cc)
    target_link_libraries(test_simple_mesh simple_mesh ${UnitTest_LIBRARIES})
    add_test(simple_mesh test_simple_mesh)

endif()

