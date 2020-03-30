# -*- mode: cmake -*-

project (INFILTRATION_1D)

add_amanzi_comparison_test(infiltration_clay_sand_1d_u
                           INPUT amanzi_infiltration_1d-u.xml
			   ABSOLUTE_TOLERANCE 1e-8
			   RELATIVE_TOLERANCE 1e-8
			   OBSERVATION "aqueous pressure"
			   REFERENCE ${INFILTRATION_1D_SOURCE_DIR}/golden_output/observations_2c.out
			   NPROCS 1)

# add_amanzi_comparison_test(infiltration_clay_sand_1d_s
#                            INPUT amanzi_infiltration_clay_sand_1d-s.xml
# 			     ABSOLUTE_TOLERANCE 1e-8
# 			     OBSERVATION "aqueous pressure"
# 			     REFERENCE observations.out
# 			     NPROCS 1)

add_amanzi_comparison_test(infiltration_loam_sand_1d_u
                           INPUT amanzi_infiltration_loam_sand_1d-u.xml
			   ABSOLUTE_TOLERANCE 1e-7
			   OBSERVATION "aqueous pressure"
			   REFERENCE ${INFILTRATION_1D_SOURCE_DIR}/golden_output/observations_2a.out
			   NPROCS 1)

# add_amanzi_comparison_test(infiltration_loam_sand_1d_s
#                            INPUT amanzi_infiltration_loam_sand_1d-s.xml
# 			     ABSOLUTE_TOLERANCE 1e-8
# 			     OBSERVATION "aqueous pressure"
# 			     REFERENCE observations.out
# 			     NPROCS 1)

add_amanzi_comparison_test(infiltration_sand_loam_1d_u
                           INPUT amanzi_infiltration_sand_loam_1d-u.xml
                           ABSOLUTE_TOLERANCE 1e-8
                           RELATIVE_TOLERANCE 1e-8
			   OBSERVATION "aqueous pressure"
			   REFERENCE ${INFILTRATION_1D_SOURCE_DIR}/golden_output/observations_2b.out
			   NPROCS 1)

# add_amanzi_comparison_test(infiltration_sand_loam_1d_s
#                            INPUT amanzi_infiltration_sand_loam_1d-s.xml
# 			     ABSOLUTE_TOLERANCE 1e-8
# 			     OBSERVATION "aqueous pressure"
# 			     REFERENCE observations.out
# 			     NPROCS 1)
