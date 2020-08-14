# -*- mode: cmake -*-

#
#  Amanzi
#   Solvers
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)

#
# Define a project name
# After this command the following varaibles are defined
#   SOLVERS_SOURCE_DIR
#   SOLVERS_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(SOLVERS)

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})
include_directories(${DBG_SOURCE_DIR})
include_directories(${WHETSTONE_SOURCE_DIR})
# -- below needed only by ResidualDebugger -- maybe this can be fixed?
include_directories(${DATA_STRUCTURES_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${OUTPUT_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${Belos_INCLUDE_DIRS})
include_directories(${NOX_INCLUDE_DIRS})
include_directories(${ASCEMIO_INCLUDE_DIRS})
include_directories(${HDF5_INCLUDE_DIRS})

#
# Library: solvers
#
file(GLOB solvers_inc_files "*.hh")
set(solvers_source_files ResidualDebugger.cc
                         PreconditionerFactory.cc
                         PreconditionerBoomerAMG.cc
                         PreconditionerML.cc
                         PreconditionerEuclid.cc
                         PreconditionerBlockILU.cc
                         SolverFactory.cc)

set(solvers_tpl_libs ${DBC_LIBRARIES}
		     ${Teuchos_LIBRARIES} ${Epetra_LIBRARIES}
                     ${Belos_LIBRARIES} ${NOX_LIBRARIES} ${Amesos2_LIBRARIES}
                     ${HYPRE_LIBRARY})
list(REVERSE solvers_tpl_libs)
list(REMOVE_DUPLICATES solvers_tpl_libs)
list(REVERSE solvers_tpl_libs)

set(solvers_link_libs data_structures
		      output
		      mesh
		      geometry
		      error_handling
		      whetstone
		      atk)

add_amanzi_library(solvers
                   SOURCE ${solvers_source_files} 
                   HEADERS ${solvers_inc_files}
		   LINK_LIBS data_structures output whetstone atk
                             ${solver_depends} ${solvers_tpl_libs} ${Epetra_LIBRARIES})


if (BUILD_TESTS)
    set(solvers_test_link_libs solvers ${solvers_link_libs})

    # Add UnitTest include directoy
    include_directories(${UnitTest_INCLUDE_DIRS})

    # Copy test directory files if an out of source build
    if (NOT (${SOLVERS_SOURCE_DIR} EQUAL ${SOLVERS_BINARY_DIR}))
        # file(GLOB DataFiles "${SOLVER_SOURCE_DIR}/test/*.xml")
        # file(COPY ${DataFiles} DESTINATION ${SOLVER_BINARY_DIR}/test/)
    endif()

    # Add the solvers directory to the include paths
    include_directories(${SOLVERS_SOURCE_DIR})

    add_amanzi_test(solver_linear_operators solver_linear_operators
                    KIND unit
                    SOURCE test/Main.cc test/solver_linear_operators.cc
                    LINK_LIBS ${solvers_test_link_libs} ${UnitTest_LIBRARIES} ${solvers_tpl_libs})

                  
    add_amanzi_test(solver_preconditioners solver_preconditioners
                    KIND unit
                    SOURCE test/Main.cc test/solver_preconditioners.cc
                    LINK_LIBS ${solvers_test_link_libs} ${UnitTest_LIBRARIES} ${solvers_tpl_libs})

    add_amanzi_test(solver_ooo_gids solver_ooo_gids
                    KIND unit
                    SOURCE test/Main.cc test/solver_ooo_gids.cc
                    LINK_LIBS ${solvers_test_link_libs} ${UnitTest_LIBRARIES} ${solvers_tpl_libs})

    add_amanzi_test(solver_linear_jfnk solver_linear_jfnk
                    KIND unit
                    SOURCE test/Main.cc test/solver_linear_jfnk.cc
                    LINK_LIBS ${solvers_test_link_libs} ${UnitTest_LIBRARIES} ${solvers_tpl_libs})
                    
    add_amanzi_test(solver_linear_nox_jfnk solver_linear_nox_jfnk
                    KIND unit
                    SOURCE test/Main.cc test/solver_linear_nox_jfnk.cc
                    LINK_LIBS ${solvers_test_link_libs} ${UnitTest_LIBRARIES} ${solvers_tpl_libs})

    add_amanzi_test(solver_nka solver_nka
                    KIND unit
                    SOURCE test/Main.cc test/solver_nka.cc
                    LINK_LIBS ${solvers_test_link_libs} ${UnitTest_LIBRARIES} ${solvers_tpl_libs})

    add_amanzi_test(solver_continuation solver_continuation
                    KIND unit
                    SOURCE test/Main.cc test/solver_continuation.cc
                    LINK_LIBS  ${solvers_test_link_libs} ${UnitTest_LIBRARIES} ${solvers_tpl_libs})

endif()
