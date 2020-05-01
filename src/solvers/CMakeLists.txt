# -*- mode: cmake -*-

#
#  Amanzi
#   Solvers
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)

project(SOLVERS)

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})
include_directories(${DBG_SOURCE_DIR})
include_directories(${WHETSTONE_SOURCE_DIR})

# -- below needed only by ResidualDebugger -- maybe this can be fixed?
include_directories(${DATA_STRUCTURES_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Tpetra_INCLUDE_DIRS})
include_directories(${Ifpack2_INCLUDE_DIRS})
include_directories(${ASCEMIO_INCLUDE_DIRS})
include_directories(${HDF5_INCLUDE_DIRS})

#
# Library: solvers
#
set(solvers_inc_files FnBaseDefs.hh
                      #ResidualDebugger.hh
                      Preconditioner.hh
                      PreconditionerIdentity.hh
                      PreconditionerDiagonal.hh
                      PreconditionerAssembled.hh
                      PreconditionerIfpack2.hh
                      PreconditionerFactory.hh
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

set(solvers_tpl_libs ${DBC_LIBRARIES}
		     ${Teuchos_LIBRARIES}
                     ${Tpetra_LIBRARIES}
                     ${Ifpack2_LIBRARIES})

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
        file(GLOB DataFiles "${SOLVERS_SOURCE_DIR}/test/*.xml")
        file(COPY ${DataFiles} DESTINATION ${SOLVERS_BINARY_DIR}/test/)
    endif()

    # Add the solvers directory to the include paths
    include_directories(${SOLVERS_SOURCE_DIR})

    add_amanzi_test(solvers_preconditioners solvers_preconditioners
                    KIND unit
                    SOURCE test/Main.cc test/solvers_preconditioners.cc
                    LINK_LIBS ${solvers_test_link_libs} ${UnitTest_LIBRARIES} ${solvers_tpl_libs})

    add_amanzi_test(solvers_linear_operators solvers_linear_operators
                    KIND unit
                    SOURCE test/Main.cc test/solvers_linear_operators.cc
                    LINK_LIBS ${solvers_test_link_libs} ${UnitTest_LIBRARIES} ${solvers_tpl_libs})

    # NOTE: this is a hypre test the checks for discontiguous GID orderings                 
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
