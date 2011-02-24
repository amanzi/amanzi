include_directories(${PROJECT_SOURCE_DIR}/utils)
include_directories(${PROJECT_SOURCE_DIR}/flow)

link_directories(${Trilinos_LIBRARY_DIRS})
add_executable(richards-nlk richards-nlk.cc)
target_link_libraries(richards-nlk moab_mesh flow cell_geometry gmv_mesh ${Trilinos_LIBRARIES})

add_executable(richards-jfnk richards-jfnk.cc)
target_link_libraries(richards-jfnk moab_mesh flow cell_geometry gmv_mesh ${Trilinos_LIBRARIES})

