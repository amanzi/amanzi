project(DUAL_POROSITY_1D)

add_amanzi_comparison_test(dual_porosity_1d
                           INPUT amanzi_dual_porosity_1d-u.xml
                           RELATIVE_TOLERANCE 1e-5
                           ABSOLUTE_TOLERANCE 1e-8
                           OBSERVATION "Tc-99 aqueous concentration"
                           REFERENCE ${DUAL_POROSITY_1D_SOURCE_DIR}/golden_output/observation.out
                           NPROCS 1)

