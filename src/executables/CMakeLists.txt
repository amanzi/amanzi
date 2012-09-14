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
include(BinaryManager)

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
include_directories(${OUTPUT_SOURCE_DIR})
include_directories(${STATE_SOURCE_DIR})
include_directories(${CHEM_SOURCE_DIR})
include_directories(${CHEMPK_SOURCE_DIR})
include_directories(${TRANSPORT_SOURCE_DIR})
include_directories(${WHETSTONE_SOURCE_DIR})
include_directories(${FLOW_SOURCE_DIR})
include_directories(${MPC_SOURCE_DIR})
include_directories(${MESH_SIMPLE_SOURCE_DIR})
include_directories(${MESH_DATA_SOURCE_DIR})
include_directories(${EXO_SOURCE_DIR})
include_directories(${TIME_INTEGRATION_SOURCE_DIR})
include_directories(${FUNCS_SOURCE_DIR})
include_directories(${MFUNCS_SOURCE_DIR})
include_directories(${PLATFORM_INTERFACE_SOURCE_DIR})


# Chemistry includes directory
get_property(CHEM_INCLUDES_DIR GLOBAL PROPERTY CHEM_INCLUDES_DIR)
include_directories(${CHEM_INCLUDES_DIR})

# The mesh frameworks
set(amanzi_mesh_libs "")
list(APPEND amanzi_mesh_libs stk_mesh exodus_reader)
if(ENABLE_STK_Mesh)
  add_definitions("-DHAVE_STK_MESH")
endif()
if(ENABLE_MOAB_Mesh)
  add_definitions("-DHAVE_MOAB_MESH")
endif()
if(ENABLE_MSTK_Mesh)
  add_definitions("-DHAVE_MSTK_MESH")
endif()

# External (TPL) include directories
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
if(ENABLE_CGNS)
    include_directories(${CGNS_INCLUDE_DIRS})
endif()    
if(ENABLE_MOAB_Mesh)
    include_directories(${MOAB_INCLUDE_DIRS})
    include_directories(${MOAB_SOURCE_DIR})
    list(APPEND amanzi_mesh_libs moab_mesh)
endif()
if(ENABLE_STK_Mesh)
    include_directories(${STK_INCLUDE_DIRS})
    include_directories(${STK_SOURCE_DIR})
endif()
include_directories(${HDF5_C_INCLUDE_DIR})
include_directories(${ASCEMIO_INCLUDE_DIR})

# Need to move this up to the root -- lpritch
add_definitions("-DUSE_MPI")

# Copy input files to build directory if out of source build
 if (NOT (EXEC_SOURCE_DIR STREQUAL EXEC_BINARY_DIR))
        execute_process(COMMAND ${CMAKE_COMMAND} -E 
          copy_directory ${EXEC_SOURCE_DIR}/inputfiles ${EXEC_BINARY_DIR}/inputfiles ) 
    endif()

add_amanzi_library(udriv
                   SOURCE amanzi_unstructured_grid_simulation_driver.cpp
                   LINK_LIBS 
                            ${Epetra_LIBRARIES} 
                            ${Teuchos_LIBRARIES}
                            mesh mesh_audit state mpc error_handling
			                ${Boost_LIBRARIES})

#
# Binary: amanzi-demo-driver
#
#add_executable(amanzi-demo-driver amanzi-demo-driver.cpp)
#set(demo_driver_tpl_libs ${Teuchos_LIBRARIES} ${Epetra_LIBRARIES} ${Boost_LIBRARIES})
#set(demo_driver_amanzi_libs time_integration mpc state simple_mesh error_handling exodus_reader ${amanzi_mesh_libs})
#set(demo_driver_pk_libs chemistry transport flow)
#target_link_libraries(amanzi-demo-driver 
#                      ${demo_driver_amanzi_libs}
#                      ${demo_driver_pk_libs}
#                      ${Trilinos_LIBRARIES}
#                      ${Boost_LIBRARIES})
#add_install_binary(amanzi-demo-driver)		  

#
# Binary: amanzi-alt-driver
#
#add_executable(amanzi-alt-driver amanzi-alt-driver.cpp)
#set(amanzi-alt-driver ${Teuchos_LIBRARIES} ${Epetra_LIBRARIES} ${Boost_LIBRARIES})
#set(alt_driver_amanzi_libs mpc state error_handling mesh_factory)
#set(alt_driver_pk_libs chemistry transport flow)
#target_link_libraries(amanzi-alt-driver
#                      ${alt_driver_amanzi_libs}
#                      ${alt_driver_pk_libs}
#                      ${Trilinos_LIBRARIES}
#                      ${Boost_LIBRARIES})
#add_install_binary(amanzi-alt-driver)

#
# Binary: amanzi-new-driver
#
#add_executable(amanzi-new-driver amanzi-new-driver.cpp)
#set(new_driver_amanzi_libs mpc state error_handling mesh_factory interface)
#set(new_driver_pk_libs chemistry transport flow)
#target_link_libraries(amanzi-new-driver
#                      ${new_driver_amanzi_libs}
#                      ${new_driver_pk_libs}
#                      ${Trilinos_LIBRARIES}
#                      ${Boost_LIBRARIES})

#
# Binary: verify_mesh
add_amanzi_executable(verify_mesh
                      SOURCE verify_mesh.cc
                      LINK_LIBS mesh_factory mesh mesh_audit ${Boost_LIBRARIES} ${Epetra_LIBRARIES})
		    		    

#
# Binary: verify_hex
# verify_hex uses deprecated boost::filesystem
# methods needs to pass version information to the
# prepoccessor
if ( ENABLE_CGNS )
    add_amanzi_executable(verify_hex
                          SOURCE verify_hex.cc
                          LINK_LIBS mesh_factory
                                    mesh mesh_audit
                                    output ${Boost_LIBRARIES} ${CGNS_LIBRARIES})
    set_target_properties(verify_hex PROPERTIES
                          COMPILE_DEFINES "${Boost_FILESYSTEM_DEFINES}")
endif()

#
# Binary: stk_mesh_test_hex
#
#if ( ENABLE_STK_Mesh AND ENABLE_CGNS )
#    add_executable(stk_mesh_test_hex stk_mesh_test_hex.cc)
#    target_link_libraries(stk_mesh_test_hex 
#                         exodus_reader mesh_data stk_mesh output mesh_audit
#                         ${Epetra_LIBRARIES} ${Teuchos_LIBRARIES} ${Boost_LIBRARIES})
#    add_install_binary(stk_mesh_test_hex)		     
#endif()

#
# Binary: verify_deform
#
if ( ENABLE_STK_Mesh )
  add_amanzi_executable(verify_deform
                        SOURCE verify_deform.cc
                        LINK_LIBS mesh_factory
                                  mesh
                                  ${amanzi_mesh_libs}
                                  output
                                  ${Epetra_LIBRARIES}
                                  ${Teuchos_LIBRARIES})
endif()                                

                                
#
#
#
add_install_shellscript(splitexo.sh)

#
# Install include files
#
add_install_include_file(amanzi_unstructured_grid_simulation_driver.hpp)
