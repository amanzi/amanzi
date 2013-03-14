# -*- mode: cmake -*-

#
#  Amanzi
#   DivGrad operators
#
# MFD matrix methods for div-grad elliptic operators

include_directories(${ATS_SOURCE_DIR}/src/)
include_directories(upwind_scheme)

add_library(divgrad matrix_mfd.cc
                    matrix_mfd_tpfa.cc
                    matrix_mfd_surf.cc
                    matrix_coupled_mfd.cc
                    upwind_scheme/upwind_cell_centered.cc
                    upwind_scheme/upwind_arithmetic_mean.cc
                    upwind_scheme/upwind_total_flux.cc
                    upwind_scheme/upwind_potential_difference.cc
                    upwind_scheme/upwind_gravity_flux.cc)

if (BUILD_TESTS)
    # Add UnitTest includes
    include_directories(${Amanzi_TPL_UnitTest_INCLUDE_DIRS})

    # Add UnitTest includes
    include_directories(${ATS_SOURCE_DIR}/src/operators/divgrad)

    add_executable(test_mfd
	  test/Main.cc
	  test/test_mfd_surf.cc)
    target_link_libraries(test_mfd
      divgrad
      amanzi_whetstone
      amanzi_geometry
      amanzi_mesh_factory
      amanzi_mstk_mesh
      amanzi_data_structures
      ${Amanzi_TPL_UnitTest_LIBRARIES} ${Amanzi_TPL_Trilinos_LIBRARIES})

endif()
