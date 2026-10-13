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
include_directories(${DATA_STRUCTURES_SOURCE_DIR})
include_directories(${OUTPUT_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Tpetra_INCLUDE_DIRS})
include_directories(${HDF5_C_INCLUDE_DIR})
include_directories(${ASCEMIO_INCLUDE_DIR})
include_directories(${NetCDF_INCLUDE_DIRS})

if (ENABLE_Silo)
  include_directories(${Silo_INCLUDE_DIRS})
endif(ENABLE_Silo)

#
# Library: output
#

set(output_tpl_list
  ${Tpetra_LIBRARIES}
  ${Teuchos_LIBRARIES}
  ${MSTK_LIBRARIES}
  ${MOAB_LIBRARIES}
  ${ASCEMIO_LIBRARIES}
  ${HDF5_LIBRARIES}
  ${NetCDF_C_LIBRARIES}
  ${Silo_LIBRARIES}
  )

# InputOutputHDF5 and InputOutputVerticallyStructuredNetCDF are both
# current.  OutputDomainSet and OutputSilo depend only on Output.hh (the
# abstract interface), not on FileHDF5/FileNetCDF directly, so they build
# against the current API too.
set(output_inc_files OutputUtils.hh
                     Output.hh
                     Input.hh
                     InputOutputHDF5.hh
                     InputOutputVerticallyStructuredNetCDF.hh
                     InputFactory.hh
                     OutputFactory.hh
                     OutputXDMF.hh
                     FileXDMF.hh
                     OutputDomainSet.hh
                     )

set(output_source_files OutputUtils.cc
                        InputOutputHDF5.cc
                        InputOutputVerticallyStructuredNetCDF.cc
                        InputFactory.cc
                        OutputFactory.cc
                        FileXDMF.cc
                        OutputXDMF.cc
                        OutputDomainSet.cc
                        )

set(output_link_libs data_structures mesh geometry atk error_handling ${output_tpl_list})


if (ENABLE_Silo)
  list(APPEND output_inc_files "OutputSilo.hh")
  list(APPEND output_source_files "OutputSilo.cc")
  add_definitions(-DENABLE_Silo)
endif(ENABLE_Silo)

add_amanzi_library(output 
                   SOURCE ${output_source_files}
                   HEADERS ${output_inc_files}
                   LINK_LIBS ${output_link_libs} ${output_tpl_list})


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
  # Add UnitTest
  include_directories(${UnitTest_INCLUDE_DIRS})

  set(output_test_link_libs output mesh_factory)
  
  if (ENABLE_MESH_MSTK)
    include_directories(${MESH_MSTK_SOURCE_DIR})
    list(APPEND output_test_link_libs mesh_mstk)
  endif()
  
  if (ENABLE_MESH_MOAB)
    include_directories(${MESH_MOAB_SOURCE_DIR})
    list(APPEND output_test_link_libs mesh_moab)
  endif()
  include_directories(${MESH_LOGICAL_SOURCE_DIR})
  include_directories(${MESH_EXTRACTED_SOURCE_DIR})
  include_directories(${MESH_FACTORY_SOURCE_DIR})
  
  # Copy test directory files if an out of source build
  if (NOT (${OUTPUT_SOURCE_DIR} EQUAL ${OUTPUT_BINARY_DIR}))
      file(GLOB DataFiles "${OUTPUT_SOURCE_DIR}/test/*.exo")
      file(GLOB DataFiles2 "${OUTPUT_SOURCE_DIR}/test/*.xmf")
      file(GLOB DataFiles3 "${OUTPUT_SOURCE_DIR}/test/*.h5")
      file(COPY ${DataFiles} ${DataFiles2} ${DataFiles3} DESTINATION ${OUTPUT_BINARY_DIR}/test/)
  endif()


  # Test: inputoutput_hdf5
  add_amanzi_test(output_inputoutput_hdf5 output_inputoutput_hdf5
                  SOURCE      test/Main.cc test/output_inputoutput_hdf5.cc
                  LINK_LIBS  ${output_test_link_libs} ${UnitTest_LIBRARIES}
                  KIND unit)
  add_amanzi_test(output_inputoutput_hdf5_np2 output_inputoutput_hdf5 NPROCS 2 KIND unit)

  # Test: inputoutput_vertically_structured_netcdf
  add_amanzi_test(output_inputoutput_vertically_structured_netcdf
                  output_inputoutput_vertically_structured_netcdf
                  SOURCE      test/Main.cc test/output_inputoutput_vertically_structured_netcdf.cc
                  LINK_LIBS  ${output_test_link_libs} ${UnitTest_LIBRARIES}
                  KIND unit)
  add_amanzi_test(output_inputoutput_vertically_structured_netcdf_np2
                  output_inputoutput_vertically_structured_netcdf NPROCS 2 KIND unit)

  if (ENABLE_MESH_MSTK OR ENABLE_MESH_MOAB)
    # Test: file_xdmf
    add_amanzi_test(output_file_xdmf output_file_xdmf
      SOURCE      test/Main.cc test/output_file_xdmf.cc
      LINK_LIBS   ${output_test_link_libs} ${UnitTest_LIBRARIES}
      KIND int)

    # Test: uses the Output interface
    add_amanzi_test(output_xdmf output_xdmf
      SOURCE         test/Main.cc test/output_xdmf.cc
      LINK_LIBS      ${output_test_link_libs} ${UnitTest_LIBRARIES}
      KIND           int)
    add_amanzi_test(output_xdmf_np2 output_xdmf NPROCS 2 KIND unit)

    if (ENABLE_Silo)
      # Test: uses the Output interface
      add_amanzi_test(output_silo output_silo
        SOURCE      test/Main.cc test/output_silo.cc
        LINK_LIBS   ${output_test_link_libs} ${UnitTest_LIBRARIES}
        KIND        int)
      add_amanzi_test(output_silo_np2 output_silo NPROCS 2 KIND unit)
    endif()

  endif()
endif(BUILD_TESTS)
