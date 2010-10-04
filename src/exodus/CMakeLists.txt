include_directories("${PROJET_SOURCE_DIR}/utils")
include_directories("${PROJECT_SOURCE_DIR}/dbc")
include_directories("${PROJECT_SOURCE_DIR}/mesh_data")

include_directories("${EXODUSII_INCLUDES}")
include_directories("${NETCDF_INCLUDES}")

add_library(exodus Exodus_readers.cc Exodus_records.cc Exodus_file.cc)
target_link_libraries(exodus dbc mesh_data ${EXODUS_LIBS_ALL})

if (BUILD_TESTS)
  add_executable(test_exodus test/Main.cc test/test_Reader.cc)
  add_test(exodus test_exodus)
  execute_process(COMMAND ${CMAKE_COMMAND} -E 
    create_symlink ${CMAKE_SOURCE_DIR}/exodus/test/test_files ${CMAKE_BINARY_DIR}/exodus/test_files) 
  target_link_libraries(test_exodus exodus unittest++)
endif()


