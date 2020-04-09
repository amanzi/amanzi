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

# compile options
set(HAVE_TRILINOS_PRECONDITIONERS FALSE)
set(HAVE_TRILINOS_SOLVERS FALSE)
set(HAVE_HYPRE_PRECONDITIONERS FALSE)

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})
include_directories(${DBG_SOURCE_DIR})
include_directories(${WHETSTONE_SOURCE_DIR})

# -- below needed only by ResidualDebugger -- maybe this can be fixed?
include_directories(${DATA_STRUCTURES_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
#include_directories(${OUTPUT_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Tpetra_INCLUDE_DIRS})
include_directories(${ASCEMIO_INCLUDE_DIRS})
include_directories(${HDF5_INCLUDE_DIRS})

if(HAVE_TRILINOS_SOLVERS)
  include_directories(${Belos_INCLUDE_DIRS})
  include_directories(${NOX_INCLUDE_DIRS})
endif()

if(HAVE_TRILINOS_PRECONDITIONERS)
  include_directories(${Ifpack_INCLUDE_DIRS})
endif()

if(HAVE_HYPRE_PRECONDITIONERS)
  include_directories(${Ifpack2_INCLUDE_DIRS})
  include_directories(${HYPRE_INCLUDE_DIRS})
endif()

#
# Library: solvers
#
set(solvers_inc_files FnBaseDefs.hh
                      #ResidualDebugger.hh
                      Preconditioner.hh
                      PreconditionerFactory.hh
                      PreconditionerDiagonal.hh
                      NKA_Base.hh
                      LinearOperator.hh
                      LinearOperatorPCG.hh
                      LinearOperatorGMRES.hh
                      LinearOperatorNKA.hh
                      Solver.hh
                      SolverNKA.hh
                      SolverFactory.hh
                      )
set(solvers_source_files SolverFactory.cc)#ResidualDebugger.cc)

if(HAVE_HYPRE_PRECONDITIONERS)
if(ENABLE_HYPRE)
  list(APPEND solvers_inc_files PreconditionerBoomerAMG.hh
                                PreconditionerEuclid.hh)

  list(APPEND solvers_source_files PreconditionerBoomerAMG.cc
                                   PreconditionerEuclid.cc)
endif()
endif()

if(HAVE_TRILINOS_PRECONDITIONERS)
  list(APPEND solvers_inc_files PreconditionerML.hh
                                PreconditionerBlockILU.hh)
  list(APPEND solvers_source_files PreconditionerML.cc
                                   PreconditionerBlockILU.cc)
endif()

set(solvers_tpl_libs ${DBC_LIBRARIES}
		     ${Teuchos_LIBRARIES}
                     ${Tpetra_LIBRARIES})

if (HAVE_HYPRE_PRECONDITIONERS)
if (ENABLE_HYPRE)
   list(APPEND solvers_tpl_libs ${HYPRE_LIBRARIES} ${Ifpack2_LIBRARIES})
endif()
endif()

if (HAVE_TRILINOS_PRECONDITIONERS)
   list(APPEND solvers_tpl_libs ${Ifpack_LIBRARIES})
endif()

list(REVERSE solvers_tpl_libs)
list(REMOVE_DUPLICATES solvers_tpl_libs)
list(REVERSE solvers_tpl_libs)

add_amanzi_library(solvers
                   SOURCE ${solvers_source_files} 
                   HEADERS ${solvers_inc_files}
		   LINK_LIBS data_structures whetstone atk ${solvers_tpl_libs})


if (BUILD_TESTS)
    set(solvers_test_link_libs data_structures
                               mesh
                               geometry
                               error_handling
                               whetstone
                               atk)

    # Add UnitTest include directoy
    include_directories(${UnitTest_INCLUDE_DIRS})

    # Copy test directory files if an out of source build
    if (NOT (${SOLVERS_SOURCE_DIR} EQUAL ${SOLVERS_BINARY_DIR}))
        # file(GLOB DataFiles "${SOLVERS_SOURCE_DIR}/test/*.xml")
        # file(COPY ${DataFiles} DESTINATION ${SOLVERS_BINARY_DIR}/test/)
    endif()

    # Add the solvers directory to the include paths
    include_directories(${SOLVERS_SOURCE_DIR})

    add_amanzi_test(solvers_preconditioners_tpetra solvers_preconditioners_tpetra
                    KIND unit
                    SOURCE test/Main.cc test/solvers_preconditioners_tpetra.cc
                    LINK_LIBS ${solvers_test_link_libs} ${UnitTest_LIBRARIES} ${solvers_tpl_libs})

    add_amanzi_test(solvers_linear_operators solvers_linear_operators
                    KIND unit
                    SOURCE test/Main.cc test/solvers_linear_operators.cc
                    LINK_LIBS ${solvers_test_link_libs} ${UnitTest_LIBRARIES} ${solvers_tpl_libs})


    # add_amanzi_test(solvers_preconditioner_noncontiguous solvers_preconditioner_noncontiguous
    #                 KIND unit
    #                 SOURCE test/Main.cc test/solvers_preconditioner_noncontiguous.cc
    #                 LINK_LIBS ${solvers_test_link_libs} ${UnitTest_LIBRARIES} ${solvers_tpl_libs})

    # add_amanzi_test(solvers_linear_jfnk solvers_linear_jfnk
    #                 KIND unit
    #                 SOURCE test/Main.cc test/solvers_linear_jfnk.cc
    #                 LINK_LIBS ${solvers_test_link_libs} ${UnitTest_LIBRARIES} ${solvers_tpl_libs})
                    
    # add_amanzi_test(solvers_linear_nox_jfnk solvers_linear_nox_jfnk
    #                 KIND unit
    #                 SOURCE test/Main.cc test/solvers_linear_nox_jfnk.cc
    #                 LINK_LIBS ${solvers_test_link_libs} ${UnitTest_LIBRARIES} ${solvers_tpl_libs})

    add_amanzi_test(solvers_nonlinear solvers_nonlinear
                    KIND unit
                    SOURCE test/Main.cc test/solvers_nonlinear.cc
                    LINK_LIBS ${solvers_test_link_libs} ${UnitTest_LIBRARIES} ${solvers_tpl_libs})

    # add_amanzi_test(solvers_continuation solvers_continuation
    #                 KIND unit
    #                 SOURCE test/Main.cc test/solvers_continuation.cc
    #                 LINK_LIBS  ${solvers_test_link_libs} ${UnitTest_LIBRARIES} ${solvers_tpl_libs})

endif()
