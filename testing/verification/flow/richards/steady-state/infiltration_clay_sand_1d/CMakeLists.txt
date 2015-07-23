# -*- mode: cmake -*-

project (INFILTRATION_CLAY_SAND_1D)

add_amanzi_comparison_test(infiltration_clay_sand_1d
                           INPUT amanzi_infiltration_clay_sand_1d-isv2.xml
			   TOLERANCE 1e-8
			   OBSERVATION "Aqueous pressure"
			   REFERENCE ${INFILTRATION_CLAY_SAND_1D_SOURCE_DIR}/golden_output/observations.out
			   NPROCS 1)


