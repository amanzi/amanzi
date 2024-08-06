# -*- mode: cmake -*-

project (HANTUSH_ANISOTROPIC_2D)

add_amanzi_comparison_test(hantush_anisotropic_2d
                           INPUT amanzi_hantush_anisotropic_2d-u.xml
                           FILES porflow4_6.exo
                           ABSOLUTE_TOLERANCE 2e-3
                           OBSERVATION drawdown
                           REFERENCE ${HANTUSH_ANISOTROPIC_2D_SOURCE_DIR}/golden_output/observation.out
                           NPROCS 1)

# add_amanzi_comparison_test(hantush_anisotropic_2d
#                            INPUT amanzi_hantush_anisotropic_2d-u.xml
#                            FILES porflow4_6.exo
#                            ABSOLUTE_TOLERANCE 2e-3
#                            FIELD hydraulic_head.cell.0
#                            REFERENCE ${HANTUSH_ANISOTROPIC_2D_SOURCE_DIR}/golden_output/plot_data.h5)
