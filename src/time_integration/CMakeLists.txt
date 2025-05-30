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
include_directories(${STATE_SOURCE_DIR})
include_directories(${STATE_SOURCE_DIR}/data)

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${NOX_INCLUDE_DIRS})
include_directories(${HDF5_C_INCLUDE_DIR})
include_directories(${ASCEMIO_INCLUDE_DIR})

#
# Library: time_integration
#
file(GLOB time_int_new_source_files "*.cc")
file(GLOB time_int_new_inc_files "*.hh")
set(time_int_new_tpl_libs ${DBC_LIBRARIES} ${Teuchos_LIBRARIES} ${Epetra_LIBRARIES} )
list(REVERSE time_int_new_tpl_libs)
list(REMOVE_DUPLICATES time_int_new_tpl_libs)
list(REVERSE time_int_new_tpl_libs)

set(time_integration_link_libs solvers state data_structures output mesh geometry error_handling whetstone atk )

add_amanzi_library(time_integration
                   SOURCE ${time_int_new_source_files}
		           HEADERS ${time_int_new_inc_files}
		           LINK_LIBS ${time_integration_link_libs} ${time_int_new_tpl_libs})

if (BUILD_TESTS)

    # Add UnitTest include directoy
    include_directories(${UnitTest_INCLUDE_DIRS})

    # Copy test directory files if an out of source build
    if (NOT (${TIME_INTEGRATION_SOURCE_DIR} EQUAL ${TIME_INTEGRATION_BINARY_DIR}))
      # file(GLOB DataFiles "${TIME_INTEGRATION_SOURCE_DIR}/test/*.xml")
      # file(COPY ${DataFiles} DESTINATION ${TIME_INTEGRATION_BINARY_DIR}/test/)
    endif()

    # Add the time_integration directory to the include paths
    include_directories(${TIME_INTEGRATION_SOURCE_DIR})
    include_directories(${TIME_INTEGRATION_SOURCE_DIR}/test)
    include_directories(${SOLVERS_SOURCE_DIR})
    include_directories(${DATA_STRUCTURES_SOURCE_DIR})

    add_amanzi_test(time_integration_solhist time_integration_solhist
                    KIND unit
		    SOURCE
		          test/Main.cc
			  test/time_integration_solhist.cc
		    LINK_LIBS
		          time_integration
			  ${time_integration_link_libs}
			  ${UnitTest_LIBRARIES}	  
			  ${time_integration_tpl_libs})

    add_amanzi_test(time_integration_nka time_integration_nka
                    KIND unit
		    SOURCE
		          test/Main.cc
			  test/time_integration_nka.cc
		    LINK_LIBS
		          time_integration
			  ${time_integration_link_libs}
			  ${UnitTest_LIBRARIES}
			  ${time_int_new_tpl_libs})

    add_amanzi_test(time_integration_rk time_integration_rk
                    KIND unit
		    SOURCE
		          test/Main.cc
			  test/time_integration_rk.cc
		    LINK_LIBS
		          time_integration
			  ${time_integration_link_libs}
			  ${UnitTest_LIBRARIES}
			  ${Teuchos_LIBRARIES}
			  ${Epetra_LIBRARIES}
			  ${time_integration_tpl_libs})

endif()
