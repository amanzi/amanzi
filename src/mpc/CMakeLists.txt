include_directories(${PROJECT_SOURCE_DIR}/utils)
include_directories(${PROJECT_SOURCE_DIR}/dbc)
include_directories(${PROJECT_SOURCE_DIR}/mesh_data)
include_directories(${PROJECT_SOURCE_DIR}/stk_mesh)
include_directories(${Trilinos_INCLUDE_DIRS})

add_library(mpc Chemistry_PK.cpp Chem_MPC.cpp DataLayout.cpp State.cpp)
target_link_libraries(mpc stk_mesh dbc mesh_data ${Trilinos_LIBRARIES})
link_directories(${Trilinos_LIBRARY_DIRS})


