# -*- mode: cmake -*-

project (UNCONFINED_SEEPAGE_1D)

add_amanzi_comparison_test(unconfined_seepage_1d
                           INPUT amanzi_unconfined_seepage_1d-u.xml
                           RELATIVE_TOLERANCE 3e-2
                           ABSOLUTE_TOLERANCE 1e-8
                           OBSERVATION "hydraulic head"
                           FILES porflow4_4.exo
                           REFERENCE ${UNCONFINED_SEEPAGE_1D_SOURCE_DIR}/golden_output/observation.out
                           NPROCS 1)
