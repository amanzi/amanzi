include_directories(${PROJECT_SOURCE_DIR}/dbc)

if (BUILD_TESTS)
  add_executable(test_boundary_conditions test/Main.cc test/test_BC.cc)
  target_link_libraries(test_boundary_conditions unittest++)
endif()
