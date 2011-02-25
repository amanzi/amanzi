# -*- mode: cmake -*-
include_directories(${PROJECT_SOURCE_DIR}/error_handling)
include_directories(${PROJECT_SOURCE_DIR}/simple_mesh)
include_directories(${PROJECT_SOURCE_DIR}/mesh_data)
include_directories(${PROJECT_SOURCE_DIR}/mesh)
include_directories(${Trilinos_INCLUDE_DIRS})
if(ENABLE_CGNS)
    include_directories(${CGNS_INCLUDES})
endif()    
add_definitions("-DUSE_MPI")
link_directories(${Trilinos_LIBRARY_DIRS})

# Define the source for the output variable
# GMV 
file(GLOB output_source_files "gmv*.c" "gmv*.cc")

if(ENABLE_CGNS)
    file(GLOB cgns_source "cgns*.cc")
    list(APPEND output_source_files "${cgns_source}")
endif()

add_library(output ${output_source_files})
set(link_libs simple_mesh error_handling mesh_data)
if (ENABLE_CGNS)
    list(APPEND ${CGNS_LIBRARY} ${HDF5_LIBS_ALL})
endif()
message("======> ${CGNS_LIBRARY}")
target_link_libraries(output ${link_libs} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})


if (BUILD_TESTS)
  add_executable(test_gmv_output test/Main.cc test/test_gmv.cc)
  target_link_libraries(test_gmv_output output unittest++ ${Trilinos_TPL_LIBRARIES})
  add_test(gmv_output test_gmv_output)
  if (ENABLE_CGNS)
      add_executable(test_cgns_output test/Main.cc test/test_cgns.cc)
      target_link_libraries(test_cgns_output output unittest++ ${CGNS_LIBRARY} ${HDF5_LIBS_ALL} 
        ${Trilinos_TPL_LIBRARIES})
      add_test(cgns_output test_cgns_output)
  endif()
endif()
