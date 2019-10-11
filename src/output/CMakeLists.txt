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
include_directories(${MSTK_SOURCE_DIR})
include_directories(${DATA_STRUCTURES_SOURCE_DIR})
include_directories(${OUTPUT_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Tpetra_INCLUDE_DIRS})
include_directories(${HDF5_C_INCLUDE_DIR})
include_directories(${ASCEMIO_INCLUDE_DIR})
include_directories(${MSTK_INCLUDE_DIR})

if (ENABLE_Silo)
  include_directories(${Silo_INCLUDE_DIRS})
endif(ENABLE_Silo)    

#
# Library: output
#

set(output_tpl_list ${Tpetra_LIBRARIES} ${Teuchos_LIBRARIES} ${MSTK_LIBRARIES} ${ASCEMIO_LIBRARIES} ${HDF5_LIBRARIES})

set(output_inc_files OutputUtils.hh
                     FileHDF5.hh
                     FileHDF5_impl.hh
                     FileXDMF.hh
                     Output.hh
                     OutputXDMF.hh
                     Input.hh
                     InputOutputHDF5.hh
                     OutputFactory.hh
                     InputFactory.hh
                     )
set(output_source_files OutputUtils.cc
                        FileHDF5.cc
                        FileXDMF.cc
                        OutputXDMF.cc
                        InputOutputHDF5.cc
                        OutputFactory.cc
                        InputFactory.cc
                        )


# if (ENABLE_Silo)
#   list(APPEND output_source_files "OutputSilo.cc")
#   list(APPEND output_tpl_list ${Silo_LIBRARIES})
#   add_definitions(-DENABLE_Silo)
# endif(ENABLE_Silo)


add_amanzi_library(output 
                   SOURCE ${output_source_files}
                   HEADERS ${output_inc_files}
                   LINK_LIBS data_structures mstk_mesh mesh geometry error_handling ${output_tpl_list})


# #
# # Binaries
# #
# add_amanzi_executable(unscramble_viz
#                       SOURCE UnscrambleViz.cc
#                       LINK_LIBS ${HDF5_LIBRARIES}
#                       OUTPUT_NAME unscramble_viz
#                       OUTPUT_DIRECTORY ${OUTPUT_BINARY_DIR})

# add_amanzi_executable(unscramble_restart
#                       SOURCE UnscrambleRestart.cc
#                       LINK_LIBS ${HDF5_LIBRARIES}
#                       OUTPUT_NAME unscramble_restart
#                       OUTPUT_DIRECTORY ${OUTPUT_BINARY_DIR})


#
# Tests
#
if (BUILD_TESTS)

  #include_directories(${MESH_FACTORY_SOURCE_DIR})

  # Add UnitTest
  include_directories(${UnitTest_INCLUDE_DIRS})

  # Copy test directory files if an out of source build
  if (NOT (${OUTPUT_SOURCE_DIR} EQUAL ${OUTPUT_BINARY_DIR}))
      file(GLOB DataFiles "${OUTPUT_SOURCE_DIR}/test/*.exo")
      file(GLOB DataFiles2 "${OUTPUT_SOURCE_DIR}/test/*.xmf")
      file(GLOB DataFiles3 "${OUTPUT_SOURCE_DIR}/test/*.h5")
      file(COPY ${DataFiles} ${DataFiles2} ${DataFiles3} DESTINATION ${OUTPUT_BINARY_DIR}/test/)
  endif()

  # # Test: gmv_output
  # add_amanzi_test(output_gmv output_gmv
  #                 SOURCE     test/Main.cc test/test_gmv.cc
  #                 LINK_LIBS  output mesh_factory ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES}
  #                 KIND unit)

  # Test: file_hdf5
  add_amanzi_test(output_file_hdf5 output_file_hdf5
                  SOURCE      test/Main.cc test/output_file_hdf5.cc
                  LINK_LIBS   output ${UnitTest_LIBRARIES}
                  KIND unit)
  add_amanzi_test(output_file_hdf5_np2 output_file_hdf5 NPROCS 2 KIND unit)
  add_amanzi_test(output_file_hdf5_np5 output_file_hdf5 NPROCS 5 KIND unit)

  # Test: file_xdmf
  add_amanzi_test(output_file_xdmf output_file_xdmf
                  SOURCE      test/Main.cc test/output_file_xdmf.cc
                  LINK_LIBS   output ${UnitTest_LIBRARIES}
                  KIND int)

  # Test: output_xdmf
  add_amanzi_test(output_xdmf output_xdmf
                  SOURCE      test/Main.cc test/output_xdmf.cc
                  LINK_LIBS   output ${UnitTest_LIBRARIES}
                  KIND int)
  add_amanzi_test(output_xdmf_np2 output_xdmf NPROCS 2 KIND unit)

  # if (ENABLE_Silo)
  #   # Test: uses the Output interface
  #   add_amanzi_test(output_silo output_silo
  #                   SOURCE      test/Main.cc test/test_silo.cc
  #                   LINK_LIBS   output mesh_factory mstk_mesh ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES}
  #                   KIND        int)
  # endif()
                  
endif(BUILD_TESTS)
