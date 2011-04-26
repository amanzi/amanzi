# -*- mode: cmake -*-

#
#  Amanzi
#    MOAB Mesh frameworks
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)

#
# Define a project name
# After this command the following varaibles are defined
#   OUTPUT_SOURCE_DIR
#   OUTPUT_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(OUTPUT)

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${SIMPLE_MESH_SOURCE_DIR})
include_directories(${MESH_DATA_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
if(ENABLE_CGNS)
    include_directories(${CGNS_INCLUDE_DIRS})
endif()    

# Need to move this up to the root level -- lpritch
add_definitions("-DUSE_MPI")

#
# Library: output
#

# GMV is always built, CGNS if enabled
file(GLOB output_source_files "gmv*.c" "gmv*.cc")
set(output_tpl_list ${Epetra_LIBRARIES} ${Teuchos_LIBRARIES})

if(ENABLE_CGNS)
    file(GLOB cgns_source "cgns*.cc")
    list(APPEND output_source_files "${cgns_source}")
    list(APPEND output_tpl_list ${CGNS_LIBRARIES})
endif()

add_library(output ${output_source_files})
target_link_libraries(output error_handling simple_mesh mesh_data ${output_tpl_list} )


if (BUILD_TESTS)

    # Add UnitTest
    include_directories(${UnitTest_INCLUDE_DIRS})

    # Test: gmv_output
    add_executable(test_gmv_output test/Main.cc test/test_gmv.cc)
    target_link_libraries(test_gmv_output output ${UnitTest_LIBRARIES})
    add_amanzi_test(gmv_output test_gmv_output KIND unit)

    # Test: cgns_output
    if (ENABLE_CGNS)
        add_executable(test_cgns_output test/Main.cc test/test_cgns.cc)
        target_link_libraries(test_cgns_output output ${UnitTest_LIBRARIES})
        add_amanzi_test(cgns_output test_cgns_output KIND int)
    endif()

endif()
