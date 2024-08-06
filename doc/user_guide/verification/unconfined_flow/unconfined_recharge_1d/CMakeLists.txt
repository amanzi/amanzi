# -*- mode: cmake -*-

project (UNCONFINED_RECHARGE_1D)

add_amanzi_comparison_test(unconfined_recharge_1d
                           INPUT amanzi_unconfined_recharge_1d-u.xml
                           ABSOLUTE_TOLERANCE 1e-8
                           RELATIVE_TOLERANCE 1e-4
                           OBSERVATION "aqueous pressure"
                           REFERENCE ${UNCONFINED_RECHARGE_1D_SOURCE_DIR}/golden_output/observation.out
                           NPROCS 1)
