include_directories(${PROJECT_SOURCE_DIR}/utils)
include_directories(${PROJECT_SOURCE_DIR}/mesh)
include_directories(${PROJECT_SOURCE_DIR}/mesh_data)
include_directories(${PROJECT_SOURCE_DIR}/moab_mesh)
include_directories(${PROJECT_SOURCE_DIR}/simple_mesh)
include_directories(${PROJECT_SOURCE_DIR}/stk_mesh)
include_directories(${PROJECT_SOURCE_DIR}/error_handling)
include_directories(${PROJECT_SOURCE_DIR}/exodus)
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
add_executable(amanzi-demo-driver amanzi-demo-driver.cpp)
target_link_libraries(amanzi-demo-driver mpc simple_mesh error_handling 
  exodus chemistry transport flow ${Trilinos_LIBRARIES})

if (ENABLE_STK_Mesh)
  target_link_libraries(amanzi-demo-driver stk_mesh)
endif()

if (ENABLE_MOAB_Mesh)
  target_link_libraries(amanzi-demo-driver ${MOAB_LIBRARIES})
endif()

if (ENABLE_MSTK_Mesh)
  target_link_libraries(amanzi-demo-driver ${MSTK_LIBRARIES})
endif()

