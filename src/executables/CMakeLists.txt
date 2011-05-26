# -*- mode: cmake -*-

#
#  Amanzi
#    Amanzi executables
#    Directory contains the source code to produce stand-alone
#    binaries.
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)

#
# Define a project name
# After this command the following varaibles are defined
#   EXEC_SOURCE_DIR
#   EXEC_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(EXEC)

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${MESH_FACTORY_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${MESH_FACTORY_SOURCE_DIR})
include_directories(${OUTPUT_SOURCE_DIR})
include_directories(${STATE_SOURCE_DIR})
include_directories(${CHEM_SOURCE_DIR})
include_directories(${TRANSPORT_SOURCE_DIR})
include_directories(${FLOW_SOURCE_DIR})
include_directories(${MPC_SOURCE_DIR})
include_directories(${MESH_SIMPLE_SOURCE_DIR})
include_directories(${MESH_DATA_SOURCE_DIR})
include_directories(${EXO_SOURCE_DIR})


# Chemistry includes directory
get_property(CHEM_INCLUDES_DIR GLOBAL PROPERTY CHEM_INCLUDES_DIR)
include_directories(${CHEM_INCLUDES_DIR})

# The mesh frameworks
set(amanzi_mesh_libs "")
list(APPEND amanzi_mesh_libs stk_mesh exodus)
add_definitions("-DHAVE_STK_MESH")
add_definitions("-DHAVE_MOAB_MESH")
add_definitions("-DHAVE_MSTK_MESH")

# External (TPL) include directories
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
if(ENABLE_CGNS)
    include_directories(${CGNS_INCLUDE_DIRS})
endif()    

# Need to move this up to the root -- lpritch
add_definitions("-DUSE_MPI")

# Copy input files to build directory if out of source build
 if (NOT (EXEC_SOURCE_DIR STREQUAL EXEC_BINARY_DIR))
        execute_process(COMMAND ${CMAKE_COMMAND} -E 
          copy_directory ${EXEC_SOURCE_DIR}/inputfiles ${EXEC_BINARY_DIR}/inputfiles ) 
    endif()

#
# Binary: amanzi-demo-driver
#
add_executable(amanzi-demo-driver amanzi-demo-driver.cpp)
set(demo_driver_tpl_libs ${Teuchos_LIBRARIES} ${Epetra_LIBRARIES} ${Boost_LIBRARIES})
set(demo_driver_amanzi_libs mpc state simple_mesh error_handling exodus ${amanzi_mesh_libs})
set(demo_driver_pk_libs chemistry transport flow)
target_link_libraries(amanzi-demo-driver 
                      ${demo_driver_amanzi_libs}
                      ${demo_driver_pk_libs}
                      ${Trilinos_LIBRARIES}
                      ${Boost_LIBRARIES})

#
# Binary: verify_mesh
#
add_executable(verify_mesh verify_mesh.cc)
target_link_libraries(verify_mesh mesh_factory mesh mesh_audit
                      ${Boost_LIBRARIES})
