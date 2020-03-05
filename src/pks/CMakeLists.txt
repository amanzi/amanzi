# -*- mode: cmake -*-

#
#  Amanzi
#    PK class
#

# Define a project name
# After this command the following varaibles are defined
#   PKS_SOURCE_DIR
#   PKS_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(PKS)

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${MESH_DATA_SOURCE_DIR})
include_directories(${MESH_FACTORY_SOURCE_DIR})
include_directories(${DATA_STRUCTURES_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${WHETSTONE_SOURCE_DIR})
include_directories(${OPERATORS_SOURCE_DIR})
include_directories(${TIME_INTEGRATION_SOURCE_DIR})
include_directories(${SOLVERS_SOURCE_DIR})
include_directories(${FUNCS_SOURCE_DIR})
include_directories(${MFUNCS_SOURCE_DIR})
include_directories(${STATE_SOURCE_DIR})
include_directories(${STATE_SOURCE_DIR}/data)
include_directories(${STATE_SOURCE_DIR}/evaluator)
include_directories(${STATE_SOURCE_DIR}/io)
include_directories(${OUTPUT_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Tpetra_INCLUDE_DIRS})
include_directories(${ASCEMIO_INCLUDE_DIR})

# linking
set(pks_tpl_libs
  ${Boost_LIBRARIES}
  ${Teuchos_LIBRARIES}
  ${Epetra_LIBRARIES}
  ${ExodusII_LIBRARIES})
list(REVERSE pks_tpl_libs)
list(REMOVE_DUPLICATES pks_tpl_libs)
list(REVERSE pks_tpl_libs)



#
# Library: pks
#
set(pks_src_files PK_Default.cc
                  PK_Factory.cc
   )

#
# Install Targets
#
set(pks_inc_files 
  PK.hh
  PK_Adaptors.hh
  PK_Default.hh
  PK_MixinLeaf.hh
  PK_MixinExplicit.hh
  PK_MixinExplicitSubcycled.hh
  PK_MixinImplicit.hh
  PK_MixinImplicitSubcycled.hh
  PK_MixinPredictorCorrector.hh  
  PK_MixinMPC.hh
  PK_MixinMPCAdvanceStepWeak.hh
  PK_MixinMPCGetDtMin.hh
  PK_MixinMPCImplicit.hh
  PK_Factory.hh
)

#
# link libraries
#
set(pks_link_libs state
  data_structures
  error_handling
  atk
  ${Tpetra_LIBRARIES})

add_amanzi_library(pks
                   SOURCE ${pks_src_files}
                   HEADERS ${pks_inc_files}
		   LINK_LIBS ${pks_link_libs})

add_install_include_file(${pks_inc_files})

#
# No tests so far
#
if (BUILD_TESTS) 

  set(pks_test_link_libs
    error_handling
    atk
    mesh
    mesh_factory
    data_structures
    mesh_functions
    state
    solvers
    time_integration
    pks
    ${UnitTest_LIBRARIES}
    )

    # Add UnitTest include directoy
    include_directories(${UnitTest_INCLUDE_DIRS})

    # Copy test directory files if an out of source build
    if (NOT (${PKS_SOURCE_DIR} EQUAL ${PKS_BINARY_DIR}))
        execute_process(COMMAND ${CMAKE_COMMAND} -E
          copy_directory ${PKS_SOURCE_DIR}/test ${PKS_BINARY_DIR}/test)
    endif()

    # Add the solvers directory to the include paths
    include_directories(${PKS_SOURCE_DIR})

    add_amanzi_test(pks_increment pks_increment
                    KIND unit
                    SOURCE test/Main.cc test/pks_increment.cc
                    LINK_LIBS ${pks_test_link_libs} ${UnitTest_LIBRARIES} ${pks_tpl_libs})

    add_amanzi_test(pks_mpc_heirarchy pks_mpc_heirarchy
                    KIND unit
                    SOURCE test/Main.cc test/pks_mpc_heirarchy.cc
                    LINK_LIBS ${pks_test_link_libs} ${UnitTest_LIBRARIES} ${pks_tpl_libs})

    add_amanzi_test(pks_mpc pks_mpc
                    KIND unit
                    SOURCE test/Main.cc test/pks_mpc.cc
                    LINK_LIBS ${pks_test_link_libs} ${UnitTest_LIBRARIES} ${pks_tpl_libs})

    add_amanzi_test(pks_ode pks_ode
                    KIND unit
                    SOURCE test/Main.cc test/pks_ode.cc
                    LINK_LIBS ${pks_test_link_libs} ${UnitTest_LIBRARIES} ${pks_tpl_libs})

    add_amanzi_test(pks_pde pks_pde
                    KIND unit
                    SOURCE test/Main.cc test/pks_pde.cc
                    LINK_LIBS ${pks_test_link_libs} ${UnitTest_LIBRARIES} ${pks_tpl_libs})
                  
endif()

# if (ENABLE_Physics)
#   add_subdirectory(eos)

#   add_subdirectory(flow)
# endif()

# # chemistry isn't REALLY physics...
# add_subdirectory(chemistry)

# # ok, transport is probably physics but oh well
# add_subdirectory(transport)

# if (ENABLE_Physics)
#   add_subdirectory(energy)
#   add_subdirectory(dummy)
# endif()

# add_subdirectory(mpc_pk)

