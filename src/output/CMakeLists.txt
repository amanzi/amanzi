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
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${MESH_SIMPLE_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${HDF5_C_INCLUDE_DIR})
include_directories(${ASCEMIO_INCLUDE_DIR})
if(ENABLE_CGNS)
    include_directories(${CGNS_INCLUDE_DIRS})
endif()    
if(ENABLE_STK_Mesh)
    add_definitions("-DHAVE_STK_MESH")
    include_directories(${STK_INCLUDE_DIRS})
    include_directories(${STK_SOURCE_DIR})
endif()


# Need to move this up to the root level -- lpritch
add_definitions("-DUSE_MPI")

#
# Library: output
#

# GMV is always built, CGNS if enabled
file(GLOB output_source_files "gmv*.c" "gmv*.cc")
set(output_tpl_list ${Epetra_LIBRARIES} ${Teuchos_LIBRARIES})

file(GLOB hdf5_source "hdf5_mesh*.cc")
list(APPEND output_source_files "${hdf5_source}")
list(APPEND output_tpl_list ${HDF5_LIBRARIES})

file(GLOB hdf5mpi_source "hdf5mpi_mesh*.cc")
list(APPEND output_source_files "${hdf5mpi_source}")
list(APPEND output_tpl_list ${ASCEMIO_LIBRARIES})

if(ENABLE_CGNS)
    file(GLOB cgns_source "cgns*.cc")
    list(APPEND output_source_files "${cgns_source}")
    list(APPEND output_tpl_list ${CGNS_LIBRARIES})
endif()

add_library(output ${output_source_files})
target_link_libraries(output simple_mesh mesh geometry error_handling ${output_tpl_list} )

#
# Install Targets
#
file(GLOB output_inc_files "*.h*")
add_install_include_file(${output_inc_files})
add_install_library(output)
#
# Install Targets
#
file(GLOB output_inc_files "*.h*")
add_install_include_file(${output_inc_files})


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
        target_link_libraries(test_cgns_output output ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES} z)
        add_amanzi_test(cgns_output test_cgns_output KIND int)
    endif()

    # Test: hdf5_output
    add_executable(test_hdf5_output test/Main.cc test/test_hdf5.cc)
    target_link_libraries(test_hdf5_output output ${UnitTest_LIBRARIES} z)
    add_test(hdf5_output test_hdf5_output)

    # Test: hdf5_mpi_output
    if (ENABLE_STK_Mesh)
        add_executable(test_hdf5mpi_output test/Main.cc test/test_hdf5mpi.cc)
        target_link_libraries(test_hdf5mpi_output stk_mesh output ${UnitTest_LIBRARIES} z)
        add_test(hdf5mpi_output test_hdf5mpi_output NPROCS 2 KIND int)
    endif()

endif()
