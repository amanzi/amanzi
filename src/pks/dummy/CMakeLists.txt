# -*- mode: cmake -*-

#
#  Amanzi
#   Dummy process kernel
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)
#include(LibraryManager)

#
# Define a project name
# After this command the following varaibles are defined
#   DUMMYPK_SOURCE_DIR
#   DUMMYPK_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(DUMMYPK)

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${MESH_DATA_SOURCE_DIR})
include_directories(${STATE_SOURCE_DIR})
include_directories(${TIME_INTEGRATION_SOURCE_DIR})
include_directories(${SOLVERS_SOURCE_DIR})
include_directories(${PKS_SOURCE_DIR})
include_directories(${MPC_SOURCE_DIR})
include_directories(${FLOW_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${HDF5_C_INCLUDE_DIR})
include_directories(${ASCEMIO_INCLUDE_DIR})

#
# Library: dummy
#
set(dummy_src_files Dummy_PK.cc)
#set(flow_tpl_libs ${Teuchos_LIBRARIES} ${NOX_LIBRARIES} ${Epetra_LIBRARIES})
#list(REVERSE flow_tpl_libs)
#list(REMOVE_DUPLICATES flow_tpl_libs)
#list(REVERSE flow_tpl_libs)
add_amanzi_library(dummy 
                   SOURCE ${dummy_src_files}
		   LINK_LIBS pks ${Epetra_LIBRARIES})
#
# Install Targets
#
#file(GLOB flow_inc_files "*.hh")
set(dummy_inc_files 
  Dummy_PK.hh
  )  

add_install_include_file(${dummy_inc_files})

# if (BUILD_TESTS) 
# #if (0)
#   include_directories(${MESH_FACTORY_SOURCE_DIR})

#   # Add UnitTest includes
#   include_directories(${UnitTest_INCLUDE_DIRS})

#   # Add MPC include directory. Need to remove this required path. -- lpritch
#   include_directories(${MPC_SOURCE_DIR})

#     # Copy test directory files if an out of source build
#     if (NOT (${DUMMYPK_SOURCE_DIR} EQUAL ${DUMMYPK_BINARY_DIR}) )
#         execute_process(COMMAND ${CMAKE_COMMAND} -E 
#           copy_directory ${DUMMYPK_SOURCE_DIR}/test ${DUMMYPK_BINARY_DIR}/test) 
#     endif()

#   set(amanzi_libs mpc flow state  mesh_factory operators 
#                     geometry mesh_functions functions dummy pks)

#   # Test: dummy_pk 
# #  add_amanzi_test(multiphase_simple multiphase_simple
# #                  KIND      int
# #                  SOURCE    test/Main.cc test/flow_multiphase_simple.cc
# #                  LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})


# endif()

