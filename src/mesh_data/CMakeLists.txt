include_directories("${PROJECT_SOURCE_DIR}/utils")
include_directories("${PROJECT_SOURCE_DIR}/dbc")

add_library(mesh_data 
  Data.cc Element_block.cc Element_field_types.cc Entity_kind.cc
  Element_types.cc Field_data.cc Node_set.cc
  Parameters.cc Side_set.cc)

target_link_libraries(mesh_data dbc)

if (BUILD_TESTS)
  add_executable(test_mesh_data 
    test/Main.cc test/test_Coords.cc test/test_Data.cc test/test_Fields.cc 
    test/test_Side_set.cc)
  add_test(test_mesh_data ${EXECUTABLE_OUTPUT_PATH}/test_mesh_data mesh_data)
  target_link_libraries(test_mesh_data mesh_data unittest++)
endif()
