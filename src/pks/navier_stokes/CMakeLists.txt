# -*- mode: cmake -*-

#
#  Amanzi
#   Viscous flow process kernel
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)

# Define a project name
# After this command the following varaibles are defined
#   NAVIER_STOKES_SOURCE_DIR
#   NAVIER_STOKES_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(NAVIER_STOKES)

# External (TPL) include directories
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${Teuchos_INCLUDE_DIRS})
# Amanzi include directories
include_directories(${ATK_SOURCE_DIR})
include_directories(${DATA_STRUCTURES_SOURCE_DIR})
include_directories(${FUNCS_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${MFUNCS_SOURCE_DIR})
include_directories(${OUTPUT_SOURCE_DIR})
include_directories(${OPERATORS_SOURCE_DIR})
include_directories(${PKS_SOURCE_DIR})
include_directories(${SOLVERS_SOURCE_DIR})
include_directories(${STATE_SOURCE_DIR})
include_directories(${TIME_INTEGRATION_SOURCE_DIR})
include_directories(${WHETSTONE_SOURCE_DIR})

#
# Library: navier_stokes
#
set(ns_src_files NavierStokes_PK.cc
                 NavierStokes_TI.cc
                 NavierStokesBoundaryFunction.cc
                 )
set(ns_tpl_libs ${Teuchos_LIBRARIES} ${Epetra_LIBRARIES} ${Boost_LIBRARIES} ${SEACAS_LIBRARIES})
if (Ifpack_ENABLE_HYPRE)
  list(APPEND ns_tpl_libs ${Ifpack_LIBRARIES} ${HYPRE_LIBRARY})
endif()
list(REVERSE ns_tpl_libs)
list(REMOVE_DUPLICATES ns_tpl_libs)
list(REVERSE ns_tpl_libs)
add_amanzi_library(navier_stokes SOURCE ${ns_src_files}
                   LINK_LIBS functions mesh_functions time_integration data_structures state
                             pks operators ${ns_tpl_libs})

#
# Install Targets
#
set(ns_inc_files 
  NavierStokes_PK.hh
  NavierStokesBoundaryFunction.hh)

add_install_include_file(${ns_inc_files})


if (BUILD_TESTS) 
    # Add UnitTest include directory
    include_directories(${UnitTest_INCLUDE_DIRS})
    include_directories(${MESH_FACTORY_SOURCE_DIR})

    # Copy test directory files if any out of source build
    if (NOT (${NAVIER_STOKES_SOURCE_DIR} EQUAL ${NAVIER_STOKES_BINARY_DIR}) )
        # execute_process(COMMAND ${CMAKE_COMMAND} -E 
        #   copy_directory ${NAVIER_STOKES_SOURCE_DIR}/test ${NAVIER_STOKES_BINARY_DIR}/test) 
        file(GLOB DataFiles "${NAVIER_STOKES_SOURCE_DIR}/test/*.xml"
                            "${NAVIER_STOKES_SOURCE_DIR}/test/*.exo")
        file(COPY ${DataFiles} DESTINATION ${NAVIER_STOKES_BINARY_DIR}/test/)
    endif()

    # Add the flow directory to the include paths
    include_directories(${NAVIER_STOKES_SOURCE_DIR})

    set(amanzi_libs geometry mesh mesh_factory state navier_stokes pks operators whetstone data_structures)

    # Test: transient 2D flow 
    add_amanzi_test(navier_stokes_2D navier_stokes_2D
                    KIND int
                    SOURCE test/Main.cc test/navier_stokes_2D.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})
endif()

