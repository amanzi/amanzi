include_directories(${PROJECT_SOURCE_DIR}/utils)
include_directories(${PROJECT_SOURCE_DIR}/dbc)
include_directories(${PROJECT_SOURCE_DIR}/mesh_data)
include_directories(${PROJECT_SOURCE_DIR}/stk_mesh)
include_directories(${PROJECT_SOURCE_DIR}/mpc)
include_directories(${Trilinos_INCLUDE_DIRS})

add_library(chemistry Chemistry_PK.cpp)
target_link_libraries(chemistry stk_mesh dbc mesh_data ${Trilinos_LIBRARIES})
link_directories(${Trilinos_LIBRARY_DIRS})


