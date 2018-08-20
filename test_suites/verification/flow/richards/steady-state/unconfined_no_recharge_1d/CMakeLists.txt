# -*- mode: cmake -*-

project (UNCONFINED_NO_RECHARGE_1D)

add_amanzi_comparison_test(unconfined_no_recharge_1d
                           INPUT amanzi_unconfined_no_recharge_1d-u.xml
                           ABSOLUTE_TOLERANCE 1e-7
                           RELATIVE_TOLERANCE 1e-4
                           OBSERVATION "aqueous pressure"
                           REFERENCE ${UNCONFINED_NO_RECHARGE_1D_SOURCE_DIR}/golden_output/observation.out
                           NPROCS 1)
