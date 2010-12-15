include_directories(${PROJECT_SOURCE_DIR}/utils)
include_directories(${PROJECT_SOURCE_DIR}/mesh_data)
include_directories(${PROJECT_SOURCE_DIR}/moab_mesh)
include_directories(${PROJECT_SOURCE_DIR}/simple_mesh)
include_directories(${PROJECT_SOURCE_DIR}/chemistry/includes)
include_directories(${PROJECT_SOURCE_DIR}/transport)
include_directories(${PROJECT_SOURCE_DIR}/flow)
include_directories(${PROJECT_SOURCE_DIR}/mpc)
include_directories(${PROJECT_SOURCE_DIR}/gmv)
include_directories(${PROJECT_SOURCE_DIR}/util)
include_directories(${Trilinos_INCLUDE_DIRS})
include_directories(${MPI_INCLUDES})
include_directories(${MOAB_INCLUDES})

add_definitions("-DUSE_MPI")

link_directories(${Trilinos_LIBRARY_DIRS})
add_executable(richards-nlk richards-nlk.cc)
target_link_libraries(richards-nlk moab_mesh simple_mesh flow utils ${Trilinos_LIBRARIES})

add_executable(richards-jfnk richards-jfnk.cc)
target_link_libraries(richards-jfnk moab_mesh simple_mesh flow utils ${Trilinos_LIBRARIES})

