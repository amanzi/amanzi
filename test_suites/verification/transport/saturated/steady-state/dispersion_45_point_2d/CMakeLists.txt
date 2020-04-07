# -*- mode: cmake -*-

project (DISPERSION_45_POINT_2D)

add_amanzi_comparison_test(dispersion_45_point_2d
                           INPUT amanzi_dispersion_45_point_2d-u.xml
                           FILES amanzi_dispersion_45_point_2d.exo
                           ABSOLUTE_TOLERANCE 5e-7
                           OBSERVATION "Tc-99 aqueous concentration"
                           REFERENCE ${DISPERSION_45_POINT_2D_SOURCE_DIR}/golden_output/obs5_2_1_r5.out
                           NPROCS 8)
