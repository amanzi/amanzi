# -*- mode: cmake -*-

project (LINEAR_FLUX_HEAD_1D)

add_amanzi_comparison_test(linear_flux_head_1d
                           INPUT amanzi_linear_flux_head_1d-u.xml
			   TOLERANCE 1e-6
			   OBSERVATION "aqueous pressure"
			   REFERENCE ${LINEAR_FLUX_HEAD_1D_SOURCE_DIR}/golden_output/observation.out
			   NPROCS 1)
