# -*- mode: cmake -*-

#
#  Amanzi
#    MOAB Mesh frameworks
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(AddParallelTest)

#
# Define a project name
# After this command the following varaibles are defined
#   MPC_SOURCE_DIR
#   MPC_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(MPC)

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${MESH_BASE_SOURCE_DIR})
include_directories(${MESH_DATA_SOURCE_DIR})
include_directories(${OUTPUT_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
#if(ENABLE_CGNS)
#    include_directories(${CGNS_INCLUDES})
#endif()    

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
add_library(state State.cpp)
set(state_tpl_list ${Boost_LIBRARIES} ${Teuchos_LIBRARIES} ${Epetra_LIBRARIES})
if (ENABLE_CGNS)
    list(APPEND state_tpl_list ${CGNS_LIBRARIES})
endif()
target_link_libraries(state output simple_mesh mesh_data cell_geometry ${state_tpl_libs})
set_property(GLOBAL PROPERTY STATE_SOURCE_DIR ${MPC_SOURCE_DIR})
set_property(GLOBAL PROPERTY STATE_BINARY_DIR ${MPC_BINARY_DIR})



#add_library(mpc MPC.cpp  State.cpp)
#if (ENABLE_CGNS)
#   target_link_libraries(mpc output simple_mesh mesh_data cell_geometry ${CGNS_LIBRARY} ${Trilinos_LIBRARIES} ${Boost_LIBRARIES})
#else(ENABLE_CGNS)
#    target_link_libraries(mpc output simple_mesh mesh_data cell_geometry ${Trilinos_LIBRARIES} ${Boost_LIBRARIES})
# endif()

#link_directories(${Trilinos_LIBRARY_DIRS})

#include_directories(${Boost_INCLUDE_DIR})
#include_directories(${Trilinos_INCLUDE_DIRS})
#include_directories(${MPI_INCLUDES})
#include_directories(${MOAB_INCLUDES})
#IF(ENABLE_CGNS)
#include_directories(${CGNS_INCLUDES})
#ENDIF()
#if (BUILD_TESTS AND ENABLE_MOAB_Mesh)
#  add_executable(test_mpc test/Main.cc test/test_driver.cpp)
#  add_test(mpc test_mpc)
#  if (NOT (CMAKE_SOURCE_DIR EQUAL CMAKE_BINARY_DIR))
#    execute_process(COMMAND ${CMAKE_COMMAND} -E 
#      copy_directory ${CMAKE_SOURCE_DIR}/mpc/test ${CMAKE_BINARY_DIR}/mpc/test) 
#  endif()
#  target_link_libraries(test_mpc mpc moab_mesh simple_mesh chemistry transport flow unittest++ ${MOAB_ALL_LIBS})
#endif()
#
