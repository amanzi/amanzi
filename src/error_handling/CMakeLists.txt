# -*- mode: cmake -*-
add_library(error_handling exceptions.cc dbc.cc errors.cc)

if (BUILD_TESTS)
  add_executable(test_error_handling
    test/Main.cc test/test_exceptions.cc test/test_dbc.cc test/test_errors.cc)
  add_test(error_handling test_error_handling)
  target_link_libraries(test_error_handling error_handling unittest++)
endif()



