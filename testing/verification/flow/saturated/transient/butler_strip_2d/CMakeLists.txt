# -*- mode: cmake -*-

project (BUTLER_STRIP_2D)

add_amanzi_comparison_test(butler_strip_2d
                           INPUT amanzi_butler_strip_2d-u.xml
                           TOLERANCE 1e-8
                           OBSERVATION Drawdown
                           REFERENCE ${BUTLER_STRIP_2D_SOURCE_DIR}/golden_output/observations.out
                           NPROCS 4)
