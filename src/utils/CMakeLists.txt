include_directories(${PROJECT_SOURCE_DIR}/gmv)
include_directories(${PROJECT_SOURCE_DIR}/dbc)
include_directories(${PROJECT_SOURCE_DIR}/mesh_data)
include_directories(${PROJECT_SOURCE_DIR}/simple_mesh)
include_directories(${Trilinos_INCLUDE_DIRS})

add_library(gmv_mesh gmv_mesh.cc)
target_link_libraries(gmv_mesh gmv simple_mesh dbc mesh_data ${Trilinos_LIBRARIES})
link_directories(${Trilinos_LIBRARY_DIRS})

if (BUILD_TESTS)
  add_executable(test_gmv_writer test/Main.cc test/test_gmv.cc)
  target_link_libraries(test_gmv_writer gmv_mesh gmv unittest++)
endif()

