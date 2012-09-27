# -*- mode: cmake -*-

#
#  Amanzi
#    State class
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)
include(LibraryManager)

#
# Define a project name
# After this command the following varaibles are defined
#   STATE_SOURCE_DIR
#   STATE_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(STATE)

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${MESH_DATA_SOURCE_DIR})
include_directories(${OUTPUT_SOURCE_DIR})
include_directories(${STK_SOURCE_DIR})
include_directories(${FUNCS_SOURCE_DIR})
include_directories(${DBG_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${ASCEMIO_INCLUDE_DIR})

# Need to move this up to the root -- lpritch
add_definitions("-DUSE_MPI")


# This state library will move to a new location 
# once Markus has created the data manager.
# I use global properties since all the PROJECT_NAME
# variables in this directory are tied to the MPC.
# These properties will be replaced with the appropriate
# directory project name once it is moved.
#
# Library: state
#
file(GLOB state_inc_files "*.hpp" "*.hh")
set(state_tpl_list ${Boost_LIBRARIES} ${Teuchos_LIBRARIES} ${Epetra_LIBRARIES})
add_amanzi_library(state 
                   SOURCE State.cpp Restart.cpp Vis.cpp chemistry_data.cc HEADERS ${state_inc_files}
                   LINK_LIBS atk functions output simple_mesh mesh_data ${state_tpl_libs})
set_property(GLOBAL PROPERTY STATE_SOURCE_DIR ${MPC_SOURCE_DIR})
set_property(GLOBAL PROPERTY STATE_BINARY_DIR ${MPC_BINARY_DIR})

#
# Install Target
#
add_install_include_file(${state_inc_files})


if (BUILD_TESTS)
  
    # Add UnitTest includes
    include_directories(${UnitTest_INCLUDE_DIRS})

    # Add State include directory. 
    include_directories(${STATE_SOURCE_DIR})

    # Copy test subdirectory for out of source builds
    if (NOT ("${STATE_SOURCE_DIR}" STREQUAL "${STATE_BINARY_DIR}"))
        execute_process(COMMAND ${CMAKE_COMMAND} -E 
          copy_directory ${STATE_SOURCE_DIR}/test ${STATE_BINARY_DIR}/test) 
    endif()
    
    #set(amanzi_libs geometry mesh state)
    
    if ( ENABLE_STK_Mesh )
      add_amanzi_test(state test_state
                      KIND int
                      SOURCE
	  	          test/Main.cc
	 		  test/test_restart.cpp
	 		  test/test_vis.cpp
                      LINK_LIBS state stk_mesh ${UnitTest_LIBRARIES})
    endif()

endif()
