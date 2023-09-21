# -*- mode: cmake -*-

project (UNCONFINED_LAYERED_2D)

add_amanzi_comparison_test(unconfined_layered_2d
                           INPUT amanzi_unconfined_layered_2d-u.xml
			   ABSOLUTE_TOLERANCE 2e-5
			   RELATIVE_TOLERANCE 1e-10
			   OBSERVATION "aqueous pressure"
			   REFERENCE ${UNCONFINED_LAYERED_2D_SOURCE_DIR}/golden_output/observation.out
			   NPROCS 1)
