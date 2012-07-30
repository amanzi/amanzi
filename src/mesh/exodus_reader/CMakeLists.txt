# -*- mode: cmake -*-

#
#  Amanzi
#    Exodus Readers
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)
include(BinaryManager)
include(LibraryManager)


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
file(GLOB exo_inc_files "*.hh")
add_amanzi_library(exodus_reader
                   SOURCE ${exo_source_files} HEADERS ${exo_inc_files}
                   LINK_LIBS 
                      error_handling mesh_data
                      ${Teuchos_LIBRARIES} ${Epetra_LIBRARIES} ${ExodusII_LIBRARIES} )


#
# Binary: exec_reader
#
add_amanzi_executable(exec_reader SOURCE test/exec_Reader.cc LINK_LIBS exodus_reader NO_INSTALL)

#
# Binary: exodus_read_test
#
add_amanzi_executable(exodus_read_test SOURCE test/exodus_read_test.cc LINK_LIBS exodus_reader NO_INSTALL)


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
    add_amanzi_test(exodus_reader test_exodus
                    KIND unit
                    SOURCE 
                         test/Main.cc 
                         test/test_Reader.cc
                         test/test_file_path.cc
                    LINK_LIBS 
                         exodus_reader
                         ${UnitTest_LIBRARIES}
                         ${Boost_FILESYSTEM_LIBRARY}
                         ${Boost_SYSTEM_LIBRARY})

    # Test: test_pexodus_2
    add_amanzi_test(pexodus_2proc test_pexodus_2
                    KIND int
                    NPROCS 2
                    SOURCE 
                         test/ParallelMain.cc 
                         test/test_pread_checkit.cc
                         test/test_file_path.cc
                    LINK_LIBS 
                          exodus_reader
                          ${UnitTest_LIBRARIES}
                          ${Boost_FILESYSTEM_LIBRARY}
                          ${Boost_SYSTEM_LIBRARY})

    # Test: test_pexodus_3
    add_amanzi_test(pexodus_3proc test_pexodus_3
                    KIND int
                    NPROCS 3
                    SOURCE 
                         test/ParallelMain.cc 
                         test/test_pread_checkit.cc
                         test/test_file_path.cc
                    LINK_LIBS 
                         exodus_reader
                         ${UnitTest_LIBRARIES}
                         ${Boost_FILESYSTEM_LIBRARY}
                         ${Boost_SYSTEM_LIBRARY})

    # Test: test_pexodus_4       
    add_amanzi_test(pexodus_4proc test_pexodus_4
                    KIND int
                    NPROCS 4
                    SOURCE 
                           test/ParallelMain.cc 
                           test/test_pread_checkit.cc
                           test/test_file_path.cc
                    LINK_LIBS 
                           exodus_reader
                           ${UnitTest_LIBRARIES}
                           ${Boost_FILESYSTEM_LIBRARY}
                           ${Boost_SYSTEM_LIBRARY})

endif()
