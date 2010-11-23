include_directories(${PROJECT_SOURCE_DIR}/error_handling)

if (BUILD_TESTS)
  add_executable(test_boundary_conditions test/Main.cc test/test_BC.cc)
  add_test(boundary_conditions test_boundary_conditions)
  target_link_libraries(test_boundary_conditions unittest++)
endif()
