# -*- mode: cmake -*-

#
#  Amanzi
#    Exodus Readers
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)


#
# Define a project name
# After this command the following varaibles are defined
#   EXO_SOURCE_DIR
#   EXO_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(EXO)

# Amanzi include directories
include_directories(${ATK_SOURCE_DIR})
include_directories(${DBC_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${MESH_DATA_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${ExodusII_INCLUDE_DIRS})

#
# Library: exodus_reader
#
file(GLOB exo_source_files "*.cc") 
add_library(exodus_reader ${exo_source_files})
add_install_library(exodus_reader)
target_link_libraries(exodus_reader 
                      error_handling mesh_data
                      ${Teuchos_LIBRARIES} ${Epetra_LIBRARIES} ${ExodusII_LIBRARIES})


#
# Binary: exec_reader
#
add_executable(exec_reader test/exec_Reader.cc)
target_link_libraries(exec_reader exodus_reader)

#
# Binary: exodus_read_test
#
add_executable(exodus_read_test test/exodus_read_test.cc)
target_link_libraries(exodus_read_test exodus_reader)

#
# Install Targets
#
file(GLOB exo_inc_files "*.hh")
add_install_include_file(${exo_inc_files})
#add_install_binary(exec_reader)
#add_install_binary(exodus_read_test)



if (BUILD_TESTS)

    # Add the UnitTest directory to the include paths
    include_directories(${UnitTest_INCLUDE_DIRS})

    # Set up the test_files directory? --lpritch
    set(EXODUS_FILE_DIR "${EXO_BINARY_DIR}/test_files")
    set_source_files_properties(test/test_file_path.cc PROPERTIES 
                                COMPILE_FLAGS
                                "-DEXODUS_FILE_DIR=\"\\\"${EXODUS_FILE_DIR}\\\"\" ")
                            
    execute_process(COMMAND ${CMAKE_COMMAND} -E 
    copy_directory ${EXO_SOURCE_DIR}/test/test_files ${EXODUS_FILE_DIR}) 

    #  On Mac OSX, compiler complains about duplicates from Boost
    #  These libraries are also listed in Epetra_LIBRARIES, Teuchos_LIBRARIES 
    #  Warnings go away when Boost_*_LIBRARY is removed. Does not impact
    #  test result may not have warnings on other systems. Error appears 
    #  for each test below. See ticket #333

    # Test: exodus_reader
    add_executable(test_exodus test/Main.cc test/test_Reader.cc test/test_file_path.cc)
    target_link_libraries(test_exodus exodus_reader ${UnitTest_LIBRARIES}
                          ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY}) 
    add_amanzi_test(exodus_reader test_exodus KIND unit)
    # set_tests_properties(exodus_reader PROPERTIES WILL_FAIL true) -- lpritch This was commented out in the original file


    # Test: test_pexodus_2
    add_executable(test_pexodus_2 test/ParallelMain.cc 
                   test/test_pread_2.cc test/test_pread_checkit.cc test/test_file_path.cc)
    target_link_libraries(test_pexodus_2 exodus_reader ${UnitTest_LIBRARIES} 
                        ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY})
    add_amanzi_test(pexodus_2proc test_pexodus_2 NPROCS 2 KIND int)

    # Test: test_pexodus_3
    add_executable(test_pexodus_3 test/ParallelMain.cc 
                   test/test_pread_3.cc test/test_pread_checkit.cc test/test_file_path.cc)
    target_link_libraries(test_pexodus_3 exodus_reader ${UnitTest_LIBRARIES} 
                          ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY})
    add_amanzi_test(pexodus_3proc test_pexodus_3 NPROCS 3 KIND int)

    # Test: test_pexodus_4       
    add_executable(test_pexodus_4 test/ParallelMain.cc 
                   test/test_pread_4.cc test/test_pread_checkit.cc test/test_file_path.cc)
    target_link_libraries(test_pexodus_4 exodus_reader  ${UnitTest_LIBRARIES}
                        ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY})
    add_amanzi_test(pexodus_4proc test_pexodus_4 NPROCS 4 KIND int)
     
endif()
