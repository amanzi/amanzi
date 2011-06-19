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
include_directories(${MESH_BASE_SOURCE_DIR})
include_directories(${MESH_DATA_SOURCE_DIR})

# mesh framework and ASCEM-IO
set(amanzi_mesh_libs "simple_mesh")
if(ENABLE_STK_Mesh)
    include_directories(${STK_INCLUDE_DIRS})
    include_directories(${STK_SOURCE_DIR})
    list(APPEND amanzi_mesh_libs stk_mesh exodus)
    add_definitions("-DHAVE_STK_MESH")
endif()

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
if(ENABLE_HDFIO)
    include_directories(${HDF5_INCLUDE_DIRS})
endif()
if(ENABLE_ASCEMIO)
    include_directories(${ASCEMIO_INCLUDE_DIRS})
endif()
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

if(ENABLE_HDFIO)
    file(GLOB hdf5_source "hdf5_mesh.cc")
    list(APPEND output_source_files "${hdf5_source}")
    list(APPEND output_tpl_list ${HDF5_LIBRARIES})
endif()

if(ENABLE_ASCEMIO)
    file(GLOB hdf5par_source "hdf5_mesh_par.cc")
    list(APPEND output_source_files "${hdf5par_source}")
    list(APPEND output_tpl_list ${ASCEMIO_LIBRARIES})
endif()

if(ENABLE_CGNS)
    file(GLOB cgns_source "cgns*.cc")
    list(APPEND output_source_files "${cgns_source}")
    list(APPEND output_tpl_list ${CGNS_LIBRARIES})
endif()

add_library(output ${output_source_files})
target_link_libraries(output error_handling ${amanzi_mesh_libs} mesh_data 
	              ${output_tpl_list} )

#
# Install Targets
#
file(GLOB output_inc_files "*.h*")
add_install_include_file(${output_inc_files})
add_install_library(output)


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

    # Test: hdf5_output
    if (ENABLE_HDFIO AND ENABLE_STK_Mesh)
        add_executable(test_hdf5_output test/Main.cc test/test_hdf5.cc)
        target_link_libraries(test_hdf5_output ${amanzi_mesh_libs} output ${UnitTest_LIBRARIES})
        add_test(hdf5_output test_hdf5_output KIND int)
    endif()

    # Test: hdf5_par_output
    if (ENABLE_ASCEMIO AND ENABLE_STK_Mesh)
        add_executable(test_hdf5par_output test/Main.cc test/test_hdf5_par.cc)
        target_link_libraries(test_hdf5par_output ${amanzi_mesh_libs} output ${UnitTest_LIBRARIES})
        add_test(hdf5par_output test_hdf5par_output NPROCS 1 KIND int)
    endif()

endif()
