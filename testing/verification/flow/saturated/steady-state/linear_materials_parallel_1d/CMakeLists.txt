# -*- mode: cmake -*-

project (LINEAR_MATERIALS_PARALLEL_1D)

add_amanzi_comparison_test(linear_materials_parallel_1d
                           INPUT amanzi_linear_materials_parallel_1d-u.xml
			   TOLERANCE 1e-7
			   OBSERVATION "aqueous pressure"
			   REFERENCE ${LINEAR_MATERIALS_PARALLEL_1D_SOURCE_DIR}/golden_output/observation.out
			   NPROCS 2)
