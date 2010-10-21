include_directories(${PROJECT_SOURCE_DIR}/utils)
include_directories(${PROJECT_SOURCE_DIR}/dbc)
include_directories(${PROJECT_SOURCE_DIR}/stk_mesh)
include_directories(${PROJECT_SOURCE_DIR}/exodus)
include_directories(${PROJECT_SOURCE_DIR}/mesh_data)

include_directories(${Trilinos_INCLUDE_DIRS})
link_directories(${Trilinos_LIBRARY_DIRS})

add_executable(read_mesh read_mesh.cc)
target_link_libraries(read_mesh dbc exodus mesh_data stk_mesh ${EXODUS_LIBS_ALL} ${Trilinos_LIBRARIES})

add_executable(stk_mesh_test_1 stk_mesh_test_1.cc)
target_link_libraries(stk_mesh_test_1 exodus mesh_data stk_mesh ${EXODUS_LIBS_ALL} ${Trilinos_LIBRARIES})







