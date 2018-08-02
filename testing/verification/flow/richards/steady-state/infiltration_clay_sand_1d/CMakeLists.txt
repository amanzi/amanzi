# -*- mode: cmake -*-

project (INFILTRATION_CLAY_SAND_1D)

add_amanzi_comparison_test(infiltration_clay_sand_1d_u
                           INPUT amanzi_infiltration_clay_sand_1d-u.xml
			   TOLERANCE 1e-8
			   OBSERVATION "aqueous pressure"
			   REFERENCE observations.out
			   NPROCS 1)

add_amanzi_comparison_test(infiltration_clay_sand_1d_s
                           INPUT amanzi_infiltration_clay_sand_1d-s.xml
			   TOLERANCE 1e-8
			   OBSERVATION "aqueous pressure"
			   REFERENCE observations.out
			   NPROCS 1)


