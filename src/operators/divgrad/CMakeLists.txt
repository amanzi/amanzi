# -*- mode: cmake -*-

#
#  Amanzi
#   DivGrad operators
#
# MFD matrix methods for div-grad elliptic operators

include_directories(upwind_scheme)

add_library(divgrad
    # MatrixMFD.cc
    # MatrixMFD_ScaledConstraint.cc
    # MatrixMFD_TPFA.cc
    # MatrixMFD_TPFA_ScaledConstraint.cc
    #                 MatrixMFD_Surf.cc
    #                 MatrixMFD_Surf_ScaledConstraint.cc
    #                 Matrix_TPFA.cc
    #                 Matrix_TPFA_Surf.cc
    #                 MatrixMFD_Coupled.cc
    #                 MatrixMFD_Coupled_TPFA.cc
    #                 MatrixMFD_Coupled_Surf.cc
    #                 MatrixMFD_Factory.cc
                    upwind_scheme/upwind_cell_centered.cc
                    upwind_scheme/upwind_arithmetic_mean.cc
                    upwind_scheme/UpwindFluxFactory.cc
                    upwind_scheme/upwind_flux_split_denominator.cc  
                    upwind_scheme/upwind_flux_fo_cont.cc 
                    upwind_scheme/upwind_flux_harmonic_mean.cc
                    upwind_scheme/upwind_total_flux.cc
                    upwind_scheme/upwind_potential_difference.cc
                    upwind_scheme/upwind_gravity_flux.cc)

install(TARGETS divgrad DESTINATION lib)

# if (BUILD_TESTS)
#     # Add UnitTest includes
#     include_directories(${Amanzi_TPL_UnitTest_INCLUDE_DIRS})

#     # Add UnitTest includes
#     include_directories(${ATS_SOURCE_DIR}/src/operators/divgrad)

#     # Copy test directory files if an out of source build
#     if (NOT (${ATS_SOURCE_DIR}/src/operators/divgrad EQUAL ${ATS_BINARY_DIR}/src/operators/divgrad) )
#         execute_process(COMMAND ${CMAKE_COMMAND} -E 
#           copy_directory ${ATS_SOURCE_DIR}/src/operators/divgrad/test ${ATS_BINARY_DIR}/src/operators/divgrad/test) 
#     endif()


#     add_executable(test_matrix_mfd test/Main.cc test/test_matrix_mfd.cc)
#     target_link_libraries(test_matrix_mfd
#       divgrad amanzi_atk amanzi_whetstone amanzi_geometry
#       amanzi_mesh_factory amanzi_mstk_mesh
#       amanzi_data_structures amanzi_solvers
#       ${Amanzi_TPL_UnitTest_LIBRARIES}
#       ${Amanzi_TPL_Trilinos_LIBRARIES})

#     add_executable(test_matrix_mfd_scaledconstraint 
#       test/Main.cc test/test_matrix_mfd_scaledconstraint.cc)
#     target_link_libraries(test_matrix_mfd_scaledconstraint
#       divgrad amanzi_atk amanzi_whetstone amanzi_geometry
#       amanzi_mesh_factory amanzi_mstk_mesh
#       amanzi_data_structures amanzi_solvers
#       ${Amanzi_TPL_UnitTest_LIBRARIES}
#       ${Amanzi_TPL_Trilinos_LIBRARIES})

#     add_executable(test_matrix_mfd_tpfa test/Main.cc test/test_matrix_mfd_tpfa.cc)
#     target_link_libraries(test_matrix_mfd_tpfa
#       divgrad amanzi_atk amanzi_whetstone amanzi_geometry
#       amanzi_mesh_factory amanzi_mstk_mesh
#       amanzi_data_structures amanzi_solvers
#       ${Amanzi_TPL_UnitTest_LIBRARIES}
#       ${Amanzi_TPL_Trilinos_LIBRARIES})

#     add_executable(test_matrix_mfd_tpfa_scaledconstraint 
#       test/Main.cc test/test_matrix_mfd_tpfa_scaledconstraint.cc)
#     target_link_libraries(test_matrix_mfd_tpfa_scaledconstraint
#       divgrad amanzi_atk amanzi_whetstone amanzi_geometry
#       amanzi_mesh_factory amanzi_mstk_mesh
#       amanzi_data_structures amanzi_solvers
#       ${Amanzi_TPL_UnitTest_LIBRARIES}
#       ${Amanzi_TPL_Trilinos_LIBRARIES})

#     add_executable(test_matrix_mfd_surf
#       test/Main.cc test/test_matrix_mfd_surf.cc)
#     target_link_libraries(test_matrix_mfd_surf
#       divgrad amanzi_atk amanzi_whetstone amanzi_geometry
#       amanzi_mesh_factory amanzi_mstk_mesh
#       amanzi_data_structures amanzi_solvers
#       ${Amanzi_TPL_UnitTest_LIBRARIES}
#       ${Amanzi_TPL_Trilinos_LIBRARIES})

#     add_executable(test_matrix_mfd_surf_sc
#       test/Main.cc test/test_matrix_mfd_surf_sc.cc)
#     target_link_libraries(test_matrix_mfd_surf_sc
#       divgrad amanzi_atk amanzi_whetstone amanzi_geometry
#       amanzi_mesh_factory amanzi_mstk_mesh
#       amanzi_data_structures amanzi_solvers
#       ${Amanzi_TPL_UnitTest_LIBRARIES}
#       ${Amanzi_TPL_Trilinos_LIBRARIES})

#     # add_executable(test_matrix_tpfa
#     #   test/Main.cc test/test_matrix_tpfa.cc)
#     # target_link_libraries(test_matrix_tpfa
#     #   divgrad amanzi_atk amanzi_whetstone amanzi_geometry
#     #   amanzi_mesh_factory amanzi_mstk_mesh
#     #   amanzi_data_structures amanzi_solvers
#     #   ${Amanzi_TPL_UnitTest_LIBRARIES}
#     #   ${Amanzi_TPL_Trilinos_LIBRARIES})

#     # add_executable(test_matrix_tpfa_surf
#     #   test/Main.cc test/test_matrix_tpfa_surf.cc)
#     # target_link_libraries(test_matrix_tpfa_surf
#     #   divgrad amanzi_atk amanzi_whetstone amanzi_geometry
#     #   amanzi_mesh_factory amanzi_mstk_mesh
#     #   amanzi_data_structures amanzi_solvers
#     #   ${Amanzi_TPL_UnitTest_LIBRARIES}
#     #   ${Amanzi_TPL_Trilinos_LIBRARIES})


#     add_executable(test_matrix_mfd_coupled
#       test/Main.cc test/test_matrix_mfd_coupled.cc)
#     target_link_libraries(test_matrix_mfd_coupled
#       divgrad amanzi_atk amanzi_whetstone amanzi_geometry
#       amanzi_mesh_factory amanzi_mstk_mesh
#       amanzi_data_structures amanzi_solvers
#       ${Amanzi_TPL_UnitTest_LIBRARIES}
#       ${Amanzi_TPL_Trilinos_LIBRARIES})

#     # add_executable(test_matrix_mfd_coupled_surf
#     #   test/Main.cc test/test_matrix_mfd_coupled_surf.cc)
#     # target_link_libraries(test_matrix_mfd_coupled_surf
#     #   divgrad amanzi_atk amanzi_whetstone amanzi_geometry
#     #   amanzi_mesh_factory amanzi_mstk_mesh
#     #   amanzi_data_structures amanzi_solvers
#     #   ${Amanzi_TPL_UnitTest_LIBRARIES}
#     #   ${Amanzi_TPL_Trilinos_LIBRARIES})

#     # add_executable(test_mfd
#     #       test/Main.cc
# 	#   test/test_mfd_surf.cc)

#     # target_link_libraries(test_mfd
#     #   divgrad
#     #   amanzi_atk
#     #   amanzi_whetstone
#     #   amanzi_geometry
#     #   amanzi_mesh_factory
#     #   amanzi_mstk_mesh
#     #   amanzi_data_structures
#     #   amanzi_solvers
#     #   ${Amanzi_TPL_UnitTest_LIBRARIES} ${Amanzi_TPL_Trilinos_LIBRARIES})


# endif()
