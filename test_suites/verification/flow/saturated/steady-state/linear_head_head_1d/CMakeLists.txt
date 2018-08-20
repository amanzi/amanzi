# -*- mode: cmake -*-

project (LINEAR_HEAD_HEAD_1D)

add_amanzi_comparison_test(linear_head_head_1d
                           INPUT amanzi_linear_head_head_1d-u.xml
			   ABSOLUTE_TOLERANCE 1e-7
			   OBSERVATION "aqueous pressure"
			   REFERENCE ${LINEAR_HEAD_HEAD_1D_SOURCE_DIR}/golden_output/observation.out
			   NPROCS 1)
