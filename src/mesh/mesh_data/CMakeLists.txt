# -*- mode: cmake -*-

#
#  Amanzi
#    Mesh Data nodes, coordinates, etc.
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)

#
# Define a project name
# After this command the following varaibles are defined
#   MESH_DATA_SOURCE_DIR
#   MESH_DATA_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(MESH_DATA)

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})

#
# Library: mesh_data
#
file(GLOB mesh_data_source_files "*.cc")
add_library(mesh_data ${mesh_data_source_files})
target_link_libraries(mesh_data error_handling ${Boost_LIBRARIES} ${Epetra_LIBRARIES})

#
# Install Targets
#
file(GLOB mesh_data_inc_files "*.hh")
add_install_include_file(${mesh_data_inc_files})
add_install_library(mesh_data)

if (BUILD_TESTS)

   # Add UnitTest include directories
   include_directories(${UnitTest_INCLUDE_DIRS})

   # Tests also use Teuchos 
   include_directories(${Teuchos_INCLUDE_DIRS})

   # Test: mesh_data
   add_executable(test_mesh_data 
                  test/Main.cc test/test_Coords.cc test/test_Data.cc test/test_Fields.cc 
                  test/test_Side_set.cc test/test_Hex_mesh.cc)
   target_link_libraries(test_mesh_data mesh_data ${UnitTest_LIBRARIES} ${Teuchos_LIBRARIES})
   add_test(mesh_data test_mesh_data KIND unit)

   # Test: mesh_hex 
   add_executable(test_hex_mesh test/Main.cc test/test_Hex_mesh.cc)
   target_link_libraries(test_hex_mesh mesh_data ${UnitTest_LIBRARIES} ${Teuchos_LIBRARIES})
   add_amanzi_test(mesh_data_parallel test_hex_mesh NPROCS 2 KIND unit)

endif()
