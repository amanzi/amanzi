# -*- mode: cmake -*-

#
#  Amanzi
#   DivGrad operators
#
# MFD matrix methods for div-grad elliptic operators

include_directories(${ATS_SOURCE_DIR}/src/)
include_directories(upwind_scheme)

add_library(divgrad MatrixMFD.cc
                    MatrixMFD_TPFA.cc
                    Matrix_TPFA.cc
                    MatrixMFD_ScaledConstraint.cc
                    MatrixMFD_TPFA_ScaledConstraint.cc
                    MatrixMFD_Surf.cc
                    MatrixMFD_Surf_ScaledConstraint.cc
                    MatrixMFD_Coupled.cc
                    MatrixMFD_Coupled_TPFA.cc
                    MatrixMFD_Coupled_Surf.cc
                    MatrixMFD_Permafrost.cc
		            MatrixMFD_Factory.cc
                    upwind_scheme/upwind_cell_centered.cc
                    upwind_scheme/upwind_arithmetic_mean.cc
                    upwind_scheme/upwind_total_flux.cc
                    upwind_scheme/upwind_potential_difference.cc
                    upwind_scheme/upwind_gravity_flux.cc)

install(TARGETS divgrad DESTINATION lib)

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
      amanzi_atk
      amanzi_whetstone
      amanzi_geometry
      amanzi_mesh_factory
      amanzi_mstk_mesh
      amanzi_data_structures
      amanzi_solvers
      ${Amanzi_TPL_UnitTest_LIBRARIES} ${Amanzi_TPL_Trilinos_LIBRARIES})

endif()
