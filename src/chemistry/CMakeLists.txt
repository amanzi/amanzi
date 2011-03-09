include_directories(includes)
include_directories(${PROJECT_SOURCE_DIR}/transport)
include_directories(${PROJECT_SOURCE_DIR}/flow)
include_directories(${PROJECT_SOURCE_DIR}/mesh)
include_directories(${PROJECT_SOURCE_DIR}/simple_mesh)
include_directories(${PROJECT_SOURCE_DIR}/mpc)

add_subdirectory(utility)
add_subdirectory(solvers)
add_subdirectory(base_chemistry)
add_subdirectory(reactions)
add_subdirectory(process_kernel)
add_subdirectory(test_drivers)

include_directories(${PROJECT_SOURCE_DIR}/mpc)
include_directories(${PROJECT_SOURCE_DIR}/transport)
include_directories(${PROJECT_SOURCE_DIR}/flow)
include_directories(${PROJECT_SOURCE_DIR}/simple_mesh)
include_directories(${PROJECT_SOURCE_DIR}/util)
include_directories(${Trilinos_INCLUDE_DIRS})

link_directories(${Trilinos_LIBRARY_DIRS})

if (BUILD_TESTS)
  add_executable(test_activity_coefficients test/Main.cc test/test_activity_coefficients.cpp)
  add_test(chemistry_activity_coefficients test_activity_coefficients)
  target_link_libraries(test_activity_coefficients chemistry unittest++)

  add_executable(test_beaker test/Main.cc test/test_beaker.cpp)
  add_test(chemistry_beaker test_beaker)
  target_link_libraries(test_beaker chemistry unittest++)

  if (NOT (CMAKE_SOURCE_DIR EQUAL CMAKE_BINARY_DIR))
    execute_process(COMMAND ${CMAKE_COMMAND} -E 
      copy_directory ${CMAKE_SOURCE_DIR}/chemistry/test ${CMAKE_BINARY_DIR}/chemistry/test) 
  endif()
endif()

