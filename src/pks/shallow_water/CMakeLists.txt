# -*- mode: cmake -*-

#
#  Amanzi
#   Shallow water process kernel
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)

# Define a project name
# After this command the following varaibles are defined
#   SHALLOW_WATER_SOURCE_DIR
#   SHALLOW_WATER_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(SHALLOW_WATER)

# External (TPL) include directories
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${Teuchos_INCLUDE_DIRS})

# Amanzi include directories
include_directories(${ATK_SOURCE_DIR})
include_directories(${DATA_STRUCTURES_SOURCE_DIR})
include_directories(${FUNCS_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${MFUNCS_SOURCE_DIR})
include_directories(${OUTPUT_SOURCE_DIR})
include_directories(${OPERATORS_SOURCE_DIR})
include_directories(${PKS_SOURCE_DIR})
include_directories(${SOLVERS_SOURCE_DIR})
include_directories(${STATE_SOURCE_DIR})
include_directories(${TIME_INTEGRATION_SOURCE_DIR})
include_directories(${WHETSTONE_SOURCE_DIR})

#
# Library: shallow_water
#
set(sw_src_files ShallowWater_PK.cc
                 ShallowWater_TI.cc
                 ShallowWaterBoundaryFunction.cc
                 DischargeEvaluator.cc
                 )
set(sw_tpl_libs ${Teuchos_LIBRARIES} ${Epetra_LIBRARIES} ${Boost_LIBRARIES} ${SEACAS_LIBRARIES})
if (Ifpack_ENABLE_HYPRE)
  list(APPEND sw_tpl_libs ${Ifpack_LIBRARIES} ${HYPRE_LIBRARY})
endif()
list(REVERSE sw_tpl_libs)
list(REMOVE_DUPLICATES sw_tpl_libs)
list(REVERSE sw_tpl_libs)
add_amanzi_library(shallow_water SOURCE ${sw_src_files}
                   LINK_LIBS functions mesh_functions time_integration data_structures state
                             pks operators ${sw_tpl_libs})

#
# Install Targets
#
set(sw_inc_files 
  ShallowWater_PK.hh
  ShallowWaterBoundaryFunction.hh)

add_install_include_file(${sw_inc_files})


if (BUILD_TESTS) 
    # Add UnitTest include directory
    include_directories(${UnitTest_INCLUDE_DIRS})
    include_directories(${MESH_FACTORY_SOURCE_DIR})

    # Copy test directory files if any out of source build
    if (NOT (${SHALLOW_WATER_SOURCE_DIR} EQUAL ${SHALLOW_WATER_BINARY_DIR}) )
        # execute_process(COMMAND ${CMAKE_COMMAND} -E 
        #   copy_directory ${SHALLOW_WATER_SOURCE_DIR}/test ${SHALLOW_WATER_BINARY_DIR}/test) 
        file(GLOB DataFiles "${SHALLOW_WATER_SOURCE_DIR}/test/*.xml"
                            "${SHALLOW_WATER_SOURCE_DIR}/test/*.exo")
        file(COPY ${DataFiles} DESTINATION ${SHALLOW_WATER_BINARY_DIR}/test/)
    endif()

    # Add the flow directory to the include paths
    include_directories(${SHALLOW_WATER_SOURCE_DIR})

    set(amanzi_libs geometry mesh mesh_factory state shallow_water pks operators whetstone data_structures)

    # Test: transient 1D flow 
    add_amanzi_test(shallow_water_1D shallow_water_1D
                    KIND int
                    NPROCS 3
                    SOURCE test/Main.cc test/shallow_water_1D.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})                
                    
    # Test: smooth 2D travelling vortex 
    add_amanzi_test(shallow_water_2D_smooth shallow_water_2D_smooth
                    KIND int
                    NPROCS 3
                    SOURCE test/Main.cc test/shallow_water_2D_smooth.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})                   
                    
    # Test: smooth 2D equilibrium
    add_amanzi_test(shallow_water_bathymetry shallow_water_bathymetry
                    KIND int
                    NPROCS 3
                    SOURCE test/Main.cc test/shallow_water_bathymetry.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})      

    # Test: lake at rest with variable bathymetry
    add_amanzi_test(shallow_water_lake_at_rest shallow_water_lake_at_rest
                    KIND int
                    NPROCS 3
                    SOURCE test/Main.cc test/shallow_water_lake_at_rest.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})
               
endif()

