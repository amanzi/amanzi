# -*- mode: cmake -*-

project (THEIS_ISOTROPIC_1D)

add_amanzi_comparison_test(theis_isotropic_1d
                           INPUT amanzi_theis_isotropic_1d-u.xml
                           ABSOLUTE_TOLERANCE 1e-8
                           OBSERVATION "aqueous pressure"
                           REFERENCE ${THEIS_ISOTROPIC_1D_SOURCE_DIR}/golden_output/observation.out
                           NPROCS 4)
