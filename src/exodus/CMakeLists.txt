include_directories("${PROJET_SOURCE_DIR}/utils")
include_directories("${PROJECT_SOURCE_DIR}/dbc")
include_directories("${PROJECT_SOURCE_DIR}/mesh_data")

include_directories("${EXODUSII_INCLUDES}")
include_directories("${NETCDF_INCLUDES}")

add_library(exodus Exodus_readers.cc Exodus_records.cc Exodus_file.cc)
target_link_libraries(exodus dbc mesh_data ${EXODUS_LIBS_ALL})

if (BUILD_TESTS)
  add_executable(test_exodus test/Main.cc test/test_Reader.cc)
  add_executable(exec_reader test/exec_Reader.cc)
  add_executable(exodus_read_test test/exodus_read_test.cc)

  add_test(exodus test_exodus)
  set_tests_properties(exodus PROPERTIES WILL_FAIL true)
  execute_process(COMMAND ${CMAKE_COMMAND} -E 
    copy_directory ${CMAKE_SOURCE_DIR}/exodus/test/test_files ${CMAKE_BINARY_DIR}/exodus/test_files) 

  target_link_libraries(test_exodus exodus unittest++ hdf5 hdf5_hl z)
  target_link_libraries(exec_reader exodus hdf5 hdf5_hl z)
  target_link_libraries(exodus_read_test exodus hdf5 hdf5_hl z)
endif()


