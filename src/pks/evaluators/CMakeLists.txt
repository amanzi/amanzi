#
#  Amanzi
#   Common evaluators
#
include(PrintVariable)
include(TestManager)

# Define a project name
project(EVALUATORS)

# External (TPL) include directories
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${Teuchos_INCLUDE_DIRS})

# Amanzi include directories
include_directories(${ATK_SOURCE_DIR})
include_directories(${DATA_STRUCTURES_SOURCE_DIR})
include_directories(${EOS_SOURCE_DIR})
include_directories(${FUNCS_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${MFUNCS_SOURCE_DIR})
include_directories(${OPERATORS_SOURCE_DIR})
include_directories(${STATE_SOURCE_DIR})

# Add library
set(evaluators_src_files 
    ApertureModelPartition.cc
    ApertureModelEvaluator.cc
    PorosityModelPartition.cc
    PorosityEvaluator.cc)
set(evaluators_tpl_libs ${Teuchos_LIBRARIES} ${Epetra_LIBRARIES} ${SEACAS_LIBRARIES})
add_amanzi_library(evaluators SOURCE ${evaluators_src_files}
                   LINK_LIBS eos functions mesh_functions data_structures state operators
                             ${evaluators_tpl_libs})

# Install targets
set(evaluators_inc_files 
    ApertureModel.hh
    ApertureModel_BartonBandis.hh
    ApertureModel_Constant.hh
    ApertureModel_ExponentialLaw.hh
    ApertureModelPartition.hh
    ApertureModelEvaluator.hh
    Porosity.hh
    PorosityModelPartition.hh)

add_install_include_file(${evaluators_inc_files})


if (BUILD_TESTS) 
    # Add UnitTest include directory
    include_directories(${UnitTest_INCLUDE_DIRS})
    include_directories(${MESH_FACTORY_SOURCE_DIR})

    # Copy test directory files if any out of source build
    if (NOT (${EVALUATORS_SOURCE_DIR} EQUAL ${EVALUATORS_BINARY_DIR}) )
        # execute_process(COMMAND ${CMAKE_COMMAND} -E 
        #   copy_directory ${EVALUATORS_SOURCE_DIR}/test ${EVALUATORS_BINARY_DIR}/test) 
        file(GLOB DataFiles "${EVALUATORS_SOURCE_DIR}/test/*.xml"
                            "${EVALUATORS_SOURCE_DIR}/test/*.exo"
                            "${EVALUATORS_SOURCE_DIR}/test/*.h5")
        file(COPY ${DataFiles} DESTINATION ${EVALUATORS_BINARY_DIR}/test/)
    endif()

    # Add the evaluators directory to the include paths
    include_directories(${EVALUATORS_SOURCE_DIR})

    set(amanzi_libs geometry eos evaluators mesh mesh_factory data_structures)

    add_amanzi_test(evaluators_aperture evaluators_aperture
                    KIND int
                    SOURCE test/Main.cc test/evaluators_aperture.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})

endif()

