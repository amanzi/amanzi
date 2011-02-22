include_directories(${PROJECT_SOURCE_DIR}/error_handling)
include_directories(${Trilinos_INCLUDE_DIRS})

add_library(cell_geometry cell_geometry.cc)

add_library(mesh_audit MeshAudit.cc)
target_link_libraries(mesh_audit error_handling cell_geometry)
link_directories(${Trilinos_LIBRARY_DIRS})

