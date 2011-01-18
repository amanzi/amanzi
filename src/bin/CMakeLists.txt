# -*- mode: cmake -*-
include_directories(${PROJECT_SOURCE_DIR}/utils)
include_directories(${PROJECT_SOURCE_DIR}/error_handling)
include_directories(${PROJECT_SOURCE_DIR}/stk_mesh)
include_directories(${PROJECT_SOURCE_DIR}/exodus)
include_directories(${PROJECT_SOURCE_DIR}/mesh_data)
include_directories(${PROJECT_SOURCE_DIR}/simple_mesh)
include_directories(${CGNS_INCLUDES})

include_directories(${Trilinos_INCLUDE_DIRS})
link_directories(${Trilinos_LIBRARY_DIRS})

if (ENABLE_STK_Mesh AND ENABLE_CGNS)

   add_executable(stk_mesh_test_hex stk_mesh_test_hex.cc)
   target_link_libraries(stk_mesh_test_hex exodus mesh_data stk_mesh
                         utils ${EXODUS_LIBS_ALL}
                         ${Trilinos_LIBRARIES} ${Boost_LIBRARIES})

endif()







