include_directories(${MPI_INCLUDES})
include_directories(${Trilinos_INCLUDE_DIRS})
include_directories(${PROJECT_SOURCE_DIR}/mpc)
include_directories(${PROJECT_SOURCE_DIR}/simple_mesh)

add_library(chemistry Chemistry_PK.cpp)
target_link_libraries(chemistry geochembase ${Trilinos_LIBRARIES})
link_directories(${Trilinos_LIBRARY_DIRS})




