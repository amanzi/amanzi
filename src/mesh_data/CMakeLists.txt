include_directories("${PROJECT_SOURCE_DIR}/utils")
include_directories("${PROJECT_SOURCE_DIR}/dbc")
include_directories(${Trilinos_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS})

add_library(mesh_data 
  Data.cc Element_block.cc Element_field_types.cc Entity_kind.cc
  Element_types.cc Field_data.cc Node_set.cc
  Parameters.cc Side_set.cc HexMeshGenerator.cc)

target_link_libraries(mesh_data dbc)

if (BUILD_TESTS)
  add_executable(test_mesh_data 
    test/Main.cc test/test_Coords.cc test/test_Data.cc test/test_Fields.cc 
    test/test_Side_set.cc test/test_Hex_mesh.cc)
  add_test(mesh_data test_mesh_data)
  target_link_libraries(test_mesh_data mesh_data unittest++
                        ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})

  add_executable(test_hex_mesh test/Main.cc test/test_Hex_mesh.cc)
  target_link_libraries(test_hex_mesh mesh_data unittest++ 
                        ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
  add_test(mesh_data_parallel 
           ${MPI_EXEC} ${MPI_EXEC_NUMPROCS_FLAG} ${MPI_EXEC_MAX_NUMPROCS} test_hex_mesh)
endif()
