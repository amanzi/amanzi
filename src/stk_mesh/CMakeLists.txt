include_directories(${PROJECT_SOURCE_DIR}/utils)
include_directories(${PROJECT_SOURCE_DIR}/dbc)
include_directories(${PROJECT_SOURCE_DIR}/mesh_data)
include_directories(${PROJECT_SOURCE_DIR}/simple_mesh)

include_directories(${Trilinos_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIR})

add_library(stk_mesh Cell_topology.cc Element_category.cc Entity_map.cc Mesh_factory.cc Mesh.cc Mesh_maps_stk.cc)
target_link_libraries(stk_mesh dbc mesh_data ${Trilinos_LIBRARIES})
link_directories(${Trilinos_LIBRARY_DIRS})

if (BUILD_TESTS)
  add_executable(test_stk_mesh test/Main.cc test/test_Mesh.cc
                 test/test_Maps.cc test/test_Hex.cc)
  add_test(stk_mesh test_stk_mesh)
  target_link_libraries(test_stk_mesh stk_mesh unittest++ ${HDF5_LIBS_ALL})
endif()

