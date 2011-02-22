# -*- mode: cmake -*-
include_directories(${PROJECT_SOURCE_DIR}/gmv)
include_directories(${PROJECT_SOURCE_DIR}/mesh)
include_directories(${PROJECT_SOURCE_DIR}/error_handling)
# include_directories(${PROJECT_SOURCE_DIR}/mesh_data)
# include_directories(${PROJECT_SOURCE_DIR}/stk_mesh)
# include_directories(${PROJECT_SOURCE_DIR}/exodus)
include_directories(${PROJECT_SOURCE_DIR}/simple_mesh)
# include_directories(${PROJECT_SOURCE_DIR}/moab_mesh)
# include_directories(${MOAB_INCLUDES})
include_directories(${Trilinos_INCLUDE_DIRS})
include_directories(${CGNS_INCLUDES})
add_definitions("-DUSE_MPI")

# GMV Mesh library
add_library(gmv_mesh gmv_mesh.cc)
target_link_libraries(gmv_mesh gmv simple_mesh error_handling mesh_data ${Trilinos_LIBRARIES})

# CGNS Library
if ( ENABLE_CGNS )
    add_library(cgns_writer cgns_mesh.cc cgns_mesh_par.cc)
    target_link_libraries(cgns_writer ${CGNS_LIBRARY} ${HDF5_LIBS_ALL} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
endif()    

# Output library
#set(out_source_files gmv_mesh.cc)
#if ( ENABLE_CGNS )
#    list(APPEND out_source_files cgns_mesh.cc cgns_mesh_par.cc)
#endif()
#add_library(output ${out_source_files})
#if (ENABLE_CGNS )
#    target_link_libraries(output ${CGNS_LIBRARY})
#endif()    

#if (NOT ENABLE_CGNS)
#    add_library(gmv_mesh gmv_mesh.cc)
#else()
#    add_library(output gmv_mesh.cc cgns_mesh.cc cgns_mesh_par.cc)
#    target_link_libraries(gmv_mesh gmv simple_mesh error_handling mesh_data ${Trilinos_LIBRARIES} ${CGNS_LIBRARY})
#endif()
link_directories(${Trilinos_LIBRARY_DIRS})

if (BUILD_TESTS)
  add_executable(test_gmv_writer test/Main.cc test/test_gmv.cc)
  target_link_libraries(test_gmv_writer gmv_mesh gmv unittest++ ${Trilinos_TPL_LIBRARIES})
  if (ENABLE_CGNS)
      add_executable(test_cgns_writer test/Main.cc test/test_cgns.cc)
      target_link_libraries(test_cgns_writer cgns_writer gmv_mesh unittest++ ${CGNS_LIBRARY} ${HDF5_LIBS_ALL} 
        ${Trilinos_TPL_LIBRARIES})
  endif()
endif()

# add_executable(verify_moab_mesh verify_moab_mesh.cc)
# target_link_libraries(verify_moab_mesh moab_mesh utils ${MOAB_ALL_LIBS})

# if (ENABLE_STK_Mesh)
#    add_executable(verify_stk_mesh verify_stk_mesh.cc)
#    target_link_libraries(verify_stk_mesh mesh_data stk_mesh exodus utils ${EXODUS_LIBS_ALL} ${Trilinos_LIBRARIES})
# endif()















