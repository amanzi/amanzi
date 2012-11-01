# -*- mode: cmake -*-
MACRO(WHOLE_ARCHIVE inlib)
  SET(${inlib} -Wl,--whole-archive ${inlib} -Wl,--no-whole-archive)
ENDMACRO(WHOLE_ARCHIVE inlib)

include_directories(${Amanzi_TPL_MSTK_INCLUDE_DIRS})
add_definitions("-DMSTK_HAVE_MPI")

add_library(pk_deform_mesh deform.cc)

if ( BUILD_TESTS )

    # Add UnitTest includes
    include_directories(${Amanzi_TPL_MSTK_INCLUDE_DIRS})
    add_definitions("-DMSTK_HAVE_MPI")

    include_directories(${Amanzi_TPL_UnitTest_INCLUDE_DIRS})
    include_directories(${ATS_SOURCE_DIR}/src/pks/flow/overland)

    # Test: miscellaneous mesh routines
    add_executable(deform_test_2D test/Main.cc test/deform_test_2D.cc )
    target_link_libraries(deform_test_2D ${pk_factory} pk_deform_mesh amanzi_mesh_functions amanzi_mesh amanzi_geometry amanzi_mesh_factory amanzi_input amanzi_output amanzi_data_structures amanzi_state ${Amanzi_TPL_UnitTest_LIBRARIES})

    add_executable(deform_test_3D test/Main.cc test/deform_test_3D.cc )
    target_link_libraries(deform_test_3D ${pk_factory} pk_deform_mesh amanzi_mesh_functions amanzi_mesh amanzi_geometry amanzi_mesh_factory amanzi_input amanzi_output amanzi_data_structures amanzi_state ${Amanzi_TPL_UnitTest_LIBRARIES})

endif()