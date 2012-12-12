# -*- mode: cmake -*-

#
#  Amanzi
#    MOAB Mesh frameworks
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)
include(BinaryManager)
include(LibraryManager)

#
# Define a project name
# After this command the following varaibles are defined
#   OUTPUT_SOURCE_DIR
#   OUTPUT_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(OUTPUT)

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${HDF5_C_INCLUDE_DIR})
include_directories(${ASCEMIO_INCLUDE_DIR})
if(ENABLE_CGNS)
    include_directories(${CGNS_INCLUDE_DIRS})
endif()    
if(ENABLE_STK_Mesh)
    add_definitions("-DHAVE_STK_MESH")
    include_directories(${STK_INCLUDE_DIRS})
    include_directories(${STK_SOURCE_DIR})
endif()

#
# Library: output
#

# GMV is always built, CGNS if enabled
file(GLOB output_source_files "gmv*.c" "gmv*.cc")
set(output_tpl_list ${Epetra_LIBRARIES} ${Teuchos_LIBRARIES})

# HDF5 serial test, HDF5 now required so always enables
file(GLOB hdf5_source "hdf5_mesh*.cc")
list(APPEND output_source_files "${hdf5_source}")
list(APPEND output_tpl_list ${HDF5_LIBRARIES})

file(GLOB hdf5mpi_source "hdf5mpi_mesh*.cc")
list(APPEND output_source_files "${hdf5mpi_source}")
list(APPEND output_tpl_list ${ASCEMIO_LIBRARIES})

if(ENABLE_CGNS)
    file(GLOB cgns_source "cgns*.cc")
    list(APPEND output_source_files "${cgns_source}")
    list(APPEND output_tpl_list ${CGNS_LIBRARIES})
endif()

file(GLOB output_inc_files "*.h*")

add_amanzi_library(output 
                   SOURCE ${output_source_files} HEADERS ${output_inc_files}
                   LINK_LIBS mesh geometry error_handling ${output_tpl_list})

foreach( lib ${output_tpl_list} )
    print_variable(lib)
    # message(STATUS "DEBUG: lib=$lib")
endforeach()

#
# Binaries
#
add_amanzi_executable(unscramble_viz
                      SOURCE unscramble_viz.cpp
                      LINK_LIBS ${HDF5_LIBRARIES}
                      OUTPUT_NAME unscramble_viz
                      OUTPUT_DIRECTORY ${OUTPUT_BINARY_DIR})

add_amanzi_executable(unscramble_restart
                      SOURCE unscramble_restart.cpp
                      LINK_LIBS ${HDF5_LIBRARIES}
                      OUTPUT_NAME unscramble_restart
                      OUTPUT_DIRECTORY ${OUTPUT_BINARY_DIR})


#
# Tests
#
if (BUILD_TESTS)

   include_directories(${MESH_FACTORY_SOURCE_DIR})

    # Add UnitTest
    include_directories(${UnitTest_INCLUDE_DIRS})

    # Test: gmv_output
    add_amanzi_test(gmv_output test_gmv_output
                    SOURCE     test/Main.cc test/test_gmv.cc
		    LINK_LIBS  output mesh_factory ${UnitTest_LIBRARIES}
		    KIND unit)

    # Test: cgns_output
    if (ENABLE_CGNS)
      add_amanzi_test(cgns_output test_cgns_output
	              SOURCE      test/Main.cc test/test_cgns.cc
		      LINK_LIBS   output mesh_factory ${UnitTest_LIBRARIES}
		      KIND int)
    endif()

    # Test: hdf5_output
    add_amanzi_test(hdf5_output test_hdf5_output
                    SOURCE      test/Main.cc test/test_hdf5.cc
		    LINK_LIBS   output mesh_factory ${UnitTest_LIBRARIES}
		    KIND int)

    # Test: hdf5_mpi_output
    if (ENABLE_STK_Mesh)
        add_amanzi_test(hdf5mpi_output test_hdf5mpi_output
                        SOURCE         test/Main.cc test/test_hdf5mpi.cc
		        LINK_LIBS      output mesh_factory ${UnitTest_LIBRARIES}
			NPROCS         2
		        KIND           int)
    endif()

endif()
