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
include_directories(${MESH_BASE_SOURCE_DIR})
include_directories(${MESH_DATA_SOURCE_DIR})
include_directories(${EXO_SOURCE_DIR})
include_directories(${OUTPUT_SOURCE_DIR})
include_directories(${STATE_SOURCE_DIR})
include_directories(${CHEM_SOURCE_DIR})
include_directories(${TRANSPORT_SOURCE_DIR})
include_directories(${FLOW_SOURCE_DIR})
include_directories(${MPC_SOURCE_DIR})

# Chemistry includes directory
get_property(CHEM_INCLUDES_DIR GLOBAL PROPERTY CHEM_INCLUDES_DIR)
include_directories(${CHEM_INCLUDES_DIR})

# The mesh frameworks
set(amanzi_mesh_libs "")
if ( ENABLE_STK_Mesh )
    include_directories(${STK_INCLUDE_DIRS})
    include_directories(${STK_SOURCE_DIR})
    list(APPEND amanzi_mesh_libs stk_mesh exodus)
    add_definitions("-DHAVE_STK_MESH")
endif()
if(ENABLE_MOAB_Mesh)
    include_directories(${MOAB_INCLUDE_DIRS})
    include_directories(${MOAB_SOURCE_DIR})
    list(APPEND amanzi_mesh_libs moab_mesh)
    add_definitions("-DHAVE_MOAB_MESH")
endif()    
if(ENABLE_MSTK_Mesh)
    include_directories(${MSTK_INCLUDE_DIRS})
    include_directories(${MSTK_SOURCE_DIR})
    list(APPEND amanzi_mesh_libs mstk_mesh)
endif()    

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
# Binary: stk_mesh_test_hex
#
if ( ENABLE_STK_Mesh AND ENABLE_CGNS )
    add_executable(stk_mesh_test_hex stk_mesh_test_hex.cc)
    target_link_libraries(stk_mesh_test_hex 
                         exodus mesh_data stk_mesh output mesh_audit
                         ${Epetra_LIBRARIES} ${Teuchos_LIBRARIES} ${Boost_LIBRARIES})
endif()

#
# Binary: verify_moab_mesh
#
if ( ENABLE_STK_Mesh )
    add_executable(verify_stk_mesh verify_stk_mesh.cc)
    target_link_libraries(verify_stk_mesh stk_mesh mesh_audit exodus)
endif()

#
# Binary: verify_moab_mesh
#
if ( ENABLE_MOAB_Mesh )
    add_executable(verify_moab_mesh verify_moab_mesh.cc)
    target_link_libraries(verify_moab_mesh moab_mesh mesh_audit)
endif()
