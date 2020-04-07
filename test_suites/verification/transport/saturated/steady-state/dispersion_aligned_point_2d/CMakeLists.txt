# -*- mode: cmake -*-

project (DISPERSION_ALIGNED_POINT_2D)

add_amanzi_comparison_test(dispersion_aligned_point_2d
                           INPUT amanzi_dispersion_aligned_point_2d-u.xml
                           FILES amanzi_dispersion_aligned_point_2d.exo
                           ABSOLUTE_TOLERANCE 1e-7
                           OBSERVATION "Tc-99 aqueous concentration"
                           REFERENCE ${DISPERSION_ALIGNED_POINT_2D_SOURCE_DIR}/golden_output/obs5_2_1_r5.out
                           NPROCS 8)
