include_directories(${PROJET_SOURCE_DIR}/utils)
include_directories(${PROJECT_SOURCE_DIR}/error_handling)
include_directories(${PROJECT_SOURCE_DIR}/mesh_data)

include_directories(${Trilinos_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIR})
include_directories(${EXODUSII_INCLUDES})
include_directories(${NETCDF_INCLUDES})

link_directories(${Trilinos_LIBRARY_DIRS})

add_library(exodus Exodus_readers.cc Exodus_records.cc Exodus_file.cc Parallel_Exodus_file.cc)
target_link_libraries(exodus error_handling mesh_data ${EXODUS_LIBRARIES})

if (BUILD_TESTS)

  set(EXODUS_FILE_DIR "${CMAKE_BINARY_DIR}/exodus/test_files")

  set_source_files_properties(test/test_file_path.cc PROPERTIES 
                              COMPILE_FLAGS
                              "-DEXODUS_FILE_DIR=\"\\\"${EXODUS_FILE_DIR}\\\"\" ")

  add_executable(test_exodus test/Main.cc test/test_Reader.cc test/test_file_path.cc)
  add_executable(exec_reader test/exec_Reader.cc)
  add_executable(exodus_read_test test/exodus_read_test.cc)

  add_test(exodus test_exodus)
  # set_tests_properties(exodus PROPERTIES WILL_FAIL true)
  execute_process(COMMAND ${CMAKE_COMMAND} -E 
    copy_directory ${CMAKE_SOURCE_DIR}/exodus/test/test_files ${EXODUS_FILE_DIR}) 

  target_link_libraries(test_exodus exodus unittest++
                        ${EXODUS_LIBRARIES}
                        ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY})
  target_link_libraries(exec_reader exodus ${EXODUS_LIBRARIES})
  target_link_libraries(exodus_read_test exodus ${EXODUS_LIBRARIES})

  add_executable(test_pexodus_2 test/ParallelMain.cc 
                 test/test_pread_2.cc test/test_pread_checkit.cc test/test_file_path.cc)
  target_link_libraries(test_pexodus_2 exodus unittest++ 
                        ${Trilinos_LIBRARIES} ${EXODUS_LIBRARIES}
                        ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${Trilinos_TPL_LIBRARIES})
  add_test(exodus_parallel_2proc 
           ${MPI_EXEC} ${MPI_EXEC_NUMPROCS_FLAG} 2 test_pexodus_2)
           
  add_executable(test_pexodus_3 test/ParallelMain.cc 
                 test/test_pread_3.cc test/test_pread_checkit.cc test/test_file_path.cc)
  target_link_libraries(test_pexodus_3 exodus unittest++ 
                        ${Trilinos_LIBRARIES} ${EXODUS_LIBRARIES}
                        ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${Trilinos_TPL_LIBRARIES})
  add_test(exodus_parallel_3proc 
           ${MPI_EXEC} ${MPI_EXEC_NUMPROCS_FLAG} 3 test_pexodus_3)
           
  add_executable(test_pexodus_4 test/ParallelMain.cc 
                 test/test_pread_4.cc test/test_pread_checkit.cc test/test_file_path.cc)
  target_link_libraries(test_pexodus_4 exodus unittest++ 
                        ${Trilinos_LIBRARIES} ${EXODUS_LIBRARIES}
                        ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${Trilinos_TPL_LIBRARIES})
  add_test(exodus_parallel_4proc 
           ${MPI_EXEC} ${MPI_EXEC_NUMPROCS_FLAG} 4 test_pexodus_4 )
           
endif()


