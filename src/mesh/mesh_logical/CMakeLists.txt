# -*- mode: cmake -*-

#
#  Amanzi
#    Mesh 
#
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)

#
# Define a project name
# After this command the following varaibles are defined
#   MESH_SOURCE_DIR
#   MESH_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(MESH_LOGICAL)

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${MESH_FACTORY_SOURCE_DIR})

# External (TPL) directories
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})

#
# Library: mesh
#
set(mesh_logical_src MeshLogical.cc
		     MeshLogicalFactory.cc
		     )

add_amanzi_library(mesh_logical SOURCE ${mesh_logical_src}
				LINK_LIBS mesh geometry ${Epetra_LIBRARIES})

#
# Install Header files
#
file(GLOB mesh_logical_inc_files "*.hh")
add_install_include_file(${mesh_logical_inc_files})


if (BUILD_TESTS)

    # Include dirs for test
    include_directories(${MESH_LOGICAL_SOURCE_DIR})
    include_directories(${GEOMETRY_SOURCE_DIR})
    include_directories(${Epetra_INCLUDE_DIRS})

    # Add UnitTest ExodusII include directories
    include_directories(${UnitTest_INCLUDE_DIRS})

    # Need to copy files for the tests 
    if (NOT (${MESH_LOGICAL_SOURCE_DIR} EQUAL ${MESH_LOGICAL_BINARY_DIR}))
      execute_process(COMMAND ${CMAKE_COMMAND} -E 
        copy_directory ${MESH_LOGICAL_SOURCE_DIR}/test ${MESH_LOGICAL_BINARY_DIR}/test) 
    endif()

    # Test link libraries
    set(test_link_libs mesh_logical mesh geometry ${UnitTest_LIBRARIES})

    # Test: init
    set(test_src_files test/Main.cc test/test_mesh_logical.cc)

    add_amanzi_test(mesh_logical test_mesh_logical
                    KIND unit
		    SOURCE ${test_src_files}
		    LINK_LIBS ${test_link_libs})

endif()
  
