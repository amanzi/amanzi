# -*- mode: cmake -*-

#
#  Amanzi
#   Operator library
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)
include(LibraryManager)

#
# Define a project name
# After this command the following varaibles are defined
#   OPERATORS_SOURCE_DIR
#   OPERATORS_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(OPERATORS)

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${MSTK_SOURCE_DIR})
include_directories(${MESH_FACTORY_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${WHETSTONE_SOURCE_DIR})
include_directories(${DATA_STRUCTURES_SOURCE_DIR})
include_directories(${SOLVERS_SOURCE_DIR})
include_directories(${SOLVERS_SOURCE_DIR})
include_directories(${OPERATORS_SOURCE_DIR})
include_directories(${OUTPUT_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${MSTK_INCLUDE_DIRS})

# Need this define. Errors from MSTK include files 
# about MPI_COMM_WORLD. --lpritch
add_definitions("-DMSTK_HAVE_MPI")

#
# Library: operator audit
#
add_amanzi_library(operator_audit SOURCE OperatorAudit.cc LINK_LIBS ${flow_tpl_libs})

#
# Library: operators
#
set(operators_inc_files 
    Operator.hh
    OperatorDefs.hh
    OperatorTypeDefs.hh
    OperatorAdvection.hh
    OperatorDiffusionFactory.hh
    OperatorDiffusion.hh
    OperatorDiffusionSurface.hh
    OperatorDiffusionWithGravity.hh
    OperatorDiffusionTPFA.hh
    OperatorSource.hh
    Reconstruction.hh
    ReconstructionCell.hh
    Upwind.hh)

add_amanzi_library(operators
                   SOURCE Operator.cc 
                   OperatorDiffusion.cc OperatorDiffusionSurface.cc 
                   OperatorDiffusionWithGravity.cc OperatorDiffusionFactory.cc
                   OperatorDiffusionTPFA.cc OperatorSource.cc 
                   OperatorAdvection.cc 
                   ReconstructionCell.cc
                   HEADERS ${operators_inc_files}
                   LINK_LIBS geometry mesh error_handling whetstone solvers data_structures ${Epetra_LIBRARIES})


if ( BUILD_TESTS )
    # Add UnitTest includes
    include_directories(${UnitTest_INCLUDE_DIRS})
    include_directories(${MESH_SIMPLE_SOURCE_DIR})

    # Copy test directory files if an out of source build
    if (NOT (${OPERATORS_SOURCE_DIR} EQUAL ${OPERATORS_BINARY_DIR}) )
        execute_process(COMMAND ${CMAKE_COMMAND} -E 
          copy_directory ${OPERATORS_SOURCE_DIR}/test ${OPERATORS_BINARY_DIR}/test) 
    endif()

    set(amanzi_libs geometry atk mesh mesh_audit mesh_factory data_structures operators whetstone output)

    # Test: matrices for diffusion
    if ( ENABLE_INDIVIDUAL_TESTS_OPERATORS OR ENABLE_INDIVIDUAL_TESTS )
        add_amanzi_test(operators_laplace_beltrami operators_laplace_beltrami
                        KIND unit
                        SOURCE test/Main.cc test/operator_laplace_beltrami.cc 
                        LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})
        add_amanzi_test(operators_parabolic_surface operators_parabolic_surface
                        KIND unit
                        SOURCE test/Main.cc test/operator_parabolic_surface.cc
                        LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})
        add_amanzi_test(operators_nonlinear_surface operators_nonlinear_surface
                        KIND unit
                        SOURCE test/Main.cc test/operator_nonlinear_surface.cc
                        LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})
    else()
        add_amanzi_test(operators_laplace_beltrami operators_laplace_beltrami
                        KIND unit
                        SOURCE test/Main.cc test/operator_laplace_beltrami.cc 
                                            test/operator_parabolic_surface.cc
                                            test/operator_nonlinear_surface.cc
                        LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})
    endif()

    add_amanzi_test(operators_diffusion operators_diffusion
                    KIND unit
                    SOURCE test/Main.cc test/operator_diffusion.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    # Test: matrices for surface
    add_amanzi_test(operators_advdiff_surface operators_advdiff_surface
                    KIND unit
                    SOURCE test/Main.cc test/operator_advdiff_surface.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    # Test: utils
    add_amanzi_test(operators_upwind operators_upwind
                    KIND unit
                    SOURCE test/Main.cc test/operator_upwind.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    # convergence analysis
    add_amanzi_test(operators_convergence operators_convergence
                    KIND unit
                    SOURCE test/Main.cc test/operator_convergence.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    add_amanzi_test(operators_convergence_2pe operators_convergence NPROCS 2 KIND unit)
endif()


