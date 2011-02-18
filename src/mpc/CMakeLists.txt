include_directories(${PROJECT_SOURCE_DIR}/utils)
include_directories(${PROJECT_SOURCE_DIR}/mesh)
include_directories(${PROJECT_SOURCE_DIR}/mesh_data)
include_directories(${PROJECT_SOURCE_DIR}/moab_mesh)
include_directories(${PROJECT_SOURCE_DIR}/simple_mesh)
include_directories(${PROJECT_SOURCE_DIR}/chemistry/includes)
include_directories(${PROJECT_SOURCE_DIR}/transport)
include_directories(${PROJECT_SOURCE_DIR}/flow)
include_directories(${PROJECT_SOURCE_DIR}/mpc)
include_directories(${PROJECT_SOURCE_DIR}/gmv)
include_directories(${PROJECT_SOURCE_DIR}/util)
include_directories(${PROJECT_SOURCE_DIR}/error_handling)
include_directories(${Boost_INCLUDE_DIR})
include_directories(${Trilinos_INCLUDE_DIRS})
include_directories(${MPI_INCLUDES})
include_directories(${MOAB_INCLUDES})
include_directories(${Boost_INCLUDES})
IF(ENABLE_CGNS)
include_directories(${CGNS_INCLUDES})
ENDIF()

add_definitions("-DUSE_MPI")




add_library(mpc MPC.cpp  State.cpp)
target_link_libraries(mpc simple_mesh gmv_mesh gmv mesh_data cell_geometry ${Trilinos_LIBRARIES} ${Boost_LIBRARIES})
link_directories(${Trilinos_LIBRARY_DIRS})

if (BUILD_TESTS AND ENABLE_MOAB_Mesh)
  add_executable(test_mpc test/Main.cc test/test_driver.cpp)
  add_test(mpc test_mpc)
  if (NOT (CMAKE_SOURCE_DIR EQUAL CMAKE_BINARY_DIR))
    execute_process(COMMAND ${CMAKE_COMMAND} -E 
      copy_directory ${CMAKE_SOURCE_DIR}/mpc/test ${CMAKE_BINARY_DIR}/mpc/test) 
  endif()
  target_link_libraries(test_mpc mpc moab_mesh simple_mesh chemistry transport flow unittest++ ${MOAB_ALL_LIBS})
endif()

