# -*- mode: cmake -*-

#
#  Amanzi
#    Time Integration class
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)

#
# Define a project name
# After this command the following varaibles are defined
#   TIME_INTEGRATION_SOURCE_DIR
#   TIME_INTEGRATION_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(TIME_INTEGRATION)

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${SOLVERS_SOURCE_DIR})
include_directories(${WHETSTONE_SOURCE_DIR})
include_directories(${DATA_STRUCTURES_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${STATE_SOURCE_DIR})
include_directories(${FUNCS_SOURCE_DIR})
include_directories(${MFUNCS_SOURCE_DIR})
include_directories(${OUTPUT_SOURCE_DIR})
include_directories(${PLATFORM_INTERFACE_SOURCE_DIR})
include_directories(${DBG_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Tpetra_INCLUDE_DIRS})
include_directories(${NOX_INCLUDE_DIRS})
include_directories(${HDF5_C_INCLUDE_DIR})
include_directories(${ASCEMIO_INCLUDE_DIR})

#
# Library: time_integration
#
set(time_integration_source_files TimestepControllerStandard.cc
                                  TimestepControllerFromFile.cc
                                  TimestepControllerSmarter.cc
                                  )
set(time_integration_inc_files TimeIntegrationDefs.hh
                               TimestepController.hh
                               TimestepControllerFixed.hh
#                               TimestepControllerAdaptive.hh
                               TimestepControllerStandard.hh
                               TimestepControllerFromFile.hh
                               TimestepControllerSmarter.hh
                               TimestepControllerFactory.hh
                               BDF_FnBase.hh
                               Explicit_TI_FnBase.hh
                               SolutionHistory.hh
                               BDF1_SolverFnBase.hh
                               BDF1_State.hh
                               BDF1_TI.hh
                               Explicit_TI_RK.hh
                               )

set(time_integration_tpl_libs ${Teuchos_LIBRARIES} ${Tpetra_LIBRARIES} )
list(REVERSE time_integration_tpl_libs)
list(REMOVE_DUPLICATES time_integration_tpl_libs)
list(REVERSE time_integration_tpl_libs)

#set(time_integration_link_libs solvers data_structures output mesh geometry error_handling whetstone atk )
set(time_integration_link_libs solvers data_structures mesh geometry error_handling atk )

add_amanzi_library(time_integration
                   SOURCE ${time_integration_source_files}
                   HEADERS ${time_integration_inc_files}
		   LINK_LIBS ${time_integration_link_libs} ${time_integration_tpl_libs})

if (BUILD_TESTS)

    # Add UnitTest include directoy
    include_directories(${UnitTest_INCLUDE_DIRS})

    # Add the time_integration directory to the include paths
    include_directories(${TIME_INTEGRATION_SOURCE_DIR})
    include_directories(${TIME_INTEGRATION_SOURCE_DIR}/test)
    include_directories(${SOLVERS_SOURCE_DIR})
    include_directories(${DATA_STRUCTURES_SOURCE_DIR})

    add_amanzi_test(time_integration_solution_history time_integration_solution_history
                    KIND unit
		    SOURCE
		          test/Main.cc
			  test/time_integration_solution_history.cc
		    LINK_LIBS
			  ${time_integration_link_libs}
			  ${UnitTest_LIBRARIES}	  
			  ${time_integration_tpl_libs})

    add_amanzi_test(time_integration_bdf1 time_integration_bdf1
                    KIND unit
        	    SOURCE
        	          test/Main.cc
        		  test/time_integration_bdf1.cc
        	    LINK_LIBS
        	          time_integration
        		  ${time_integration_link_libs}
        		  ${UnitTest_LIBRARIES}
        		  ${time_integration_tpl_libs})

    add_amanzi_test(time_integration_rk time_integration_rk
                    KIND unit
        	    SOURCE
        	          test/Main.cc
        		  test/time_integration_explicit_RK.cc
        	    LINK_LIBS
        	          time_integration
        		  ${time_integration_link_libs}
        		  ${UnitTest_LIBRARIES}
        		  ${Teuchos_LIBRARIES}
        		  ${Tpetra_LIBRARIES}
        		  ${time_integration_tpl_libs})

endif()
