# -*- mode: cmake -*-

#
#  Amanzi
#    Data Structures
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)

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
include_directories(${MESH_SOURCE_DIR})
include_directories(${DATA_STRUCTURES_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${Thyra_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})


# Library: data_structures
#
file(GLOB data_structures_src_files "*.cc")
file(GLOB data_structures_inc_files "*.hh")
set(data_structures_tpl_list ${Teuchos_LIBRARIES} ${Epetra_LIBRARIES} ${Thyra_LIBRARIES})
add_amanzi_library(data_structures 
                   SOURCE ${data_structures_src_files}
		   HEADERS ${data_structures_inc_files}
		   LINK_LIBS 
		            mesh
			    ${data_structures_tpl_list})
set_property(GLOBAL PROPERTY DATA_STRUCTURES_SOURCE_DIR ${MPC_SOURCE_DIR})
set_property(GLOBAL PROPERTY DATA_STRUCTURES_BINARY_DIR ${MPC_BINARY_DIR})

#
# Install Target
#

if (BUILD_TESTS)
    # Add UnitTest includes
    include_directories(${UnitTest_INCLUDE_DIRS})
    include_directories(${Trilinos_INCLUDE_DIRS})

    # Add Data_structures include directory.
    include_directories(${DATA_STRUCTURES_SOURCE_DIR})

    include_directories(${MESH_SOURCE_DIR})
    include_directories(${MESH_FACTORY_SOURCE_DIR})
    include_directories(${MESH_SIMPLE_SOURCE_DIR})

    # Copy test subdirectory for out of source builds
    if (NOT (DATA_STRUCTURES_SOURCE_DIR STREQUAL DATA_STRUCTURES_BINARY_DIR))
        execute_process(COMMAND ${CMAKE_COMMAND} -E
          copy_directory ${DATA_STRUCTURES_SOURCE_DIR}/test ${DATA_STRUCTURES_BINARY_DIR}/test)
    endif()

    add_amanzi_test(data_structures test_data_structures 
                    KIND int
		    SOURCE 
		          test/Main.cc 
			  test/test_composite_vector.cc
			  test/test_tree_vector.cc
			  test/test_thyra_wrappers.cc
		    LINK_LIBS
		             data_structures
			     mesh
			     geometry
			     mesh_factory
			     simple_mesh
			     ${UnitTest_LIBRARIES})

#    add_amanzi_test(test_access test_access_times
#                    KIND int
#		    SOURCE test/Main.cc test/test_cv_access.cc
#		    LINK_LIBS
#		             data_structures
#			     mesh
#			     geometry
#			     mesh_factory
#			     simple_mesh
#			     ${UnitTest_LIBRARIES})

endif()
