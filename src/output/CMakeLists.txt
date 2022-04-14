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
include_directories(${MESH_DATA_SOURCE_DIR})
include_directories(${MESH_MSTK_SOURCE_DIR})
include_directories(${OUTPUT_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${HDF5_C_INCLUDE_DIR})
include_directories(${ASCEMIO_INCLUDE_DIR})

if (ENABLE_Silo)
  include_directories(${Silo_INCLUDE_DIRS})
endif(ENABLE_Silo)    

#
# Library: output
#

# GMV is always built
file(GLOB output_source_files "gmv*.c" "GMV*.cc")
set(output_tpl_list ${Epetra_LIBRARIES} ${Teuchos_LIBRARIES})

if (ENABLE_Silo)
  list(APPEND output_source_files "OutputSilo.cc")
  list(APPEND output_tpl_list ${Silo_LIBRARIES})
  add_definitions(-DENABLE_Silo)
endif(ENABLE_Silo)

# HDF5 serial test, HDF5 now required so always enables
file(GLOB hdf5_source "HDF5Mesh*.cc")
list(APPEND output_source_files "${hdf5_source}")
list(APPEND output_tpl_list ${HDF5_LIBRARIES})

file(GLOB hdf5mpi_source "HDF5_MPI*.cc")
list(APPEND output_source_files "${hdf5mpi_source}")
list(APPEND output_source_files "OutputXDMF.cc")
list(APPEND output_tpl_list ${ASCEMIO_LIBRARIES})

file(GLOB output_inc_files "*.h*")

add_amanzi_library(output 
                   SOURCE ${output_source_files} HEADERS ${output_inc_files}
                   LINK_LIBS mesh geometry error_handling ${output_tpl_list})

foreach( lib ${output_tpl_list} )
  #print_variable(lib)
  # message(STATUS "DEBUG: lib=$lib")
endforeach()

#
# Binaries
#
add_amanzi_executable(unscramble_viz
                      SOURCE UnscrambleViz.cc
                      LINK_LIBS ${HDF5_LIBRARIES}
                      OUTPUT_NAME unscramble_viz
                      OUTPUT_DIRECTORY ${OUTPUT_BINARY_DIR})

add_amanzi_executable(unscramble_restart
                      SOURCE UnscrambleRestart.cc
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

  # Copy test directory files if an out of source build
  if (NOT (${OUTPUT_SOURCE_DIR} EQUAL ${OUTPUT_BINARY_DIR}))
      file(GLOB DataFiles "${OUTPUT_SOURCE_DIR}/test/*.exo")
      file(COPY ${DataFiles} DESTINATION ${OUTPUT_BINARY_DIR}/test/)
  endif()

    
  # Test: gmv_output
  add_amanzi_test(output_gmv output_gmv
                  SOURCE     test/Main.cc test/test_gmv.cc
                  LINK_LIBS  output mesh_factory ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES}
                  KIND unit)

  # Test: hdf5_output
  add_amanzi_test(output_hdf5 output_hdf5
                  SOURCE      test/Main.cc test/test_hdf5.cc
                  LINK_LIBS   output mesh_factory ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES}
                  KIND int)



  if (ENABLE_MESH_MSTK)

    # Test: hdf5_mpi_output
    add_amanzi_test(output_hdf5mpi output_hdf5mpi
                  SOURCE         test/Main.cc test/test_hdf5mpi.cc
                  LINK_LIBS      output mesh_factory mesh_mstk ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES}
                  NPROCS         2
                  KIND           int)

    # Test: uses the Output interface
    add_amanzi_test(output_xdmf output_xdmf
                  SOURCE         test/Main.cc test/test_xdmf.cc
                  LINK_LIBS      output mesh_factory mesh_mstk ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES}
                  KIND           int)

    if (ENABLE_Silo)
      # Test: uses the Output interface
      add_amanzi_test(output_silo output_silo
                    SOURCE      test/Main.cc test/test_silo.cc
                    LINK_LIBS   output mesh_factory mesh_mstk ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES}
                    KIND        int)
    endif()

  endif()

 
                  
endif(BUILD_TESTS)
