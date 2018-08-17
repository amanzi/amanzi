# -*- mode: cmake -*-

project (INFILTRATION_LOAM_SAND_1D)

add_amanzi_comparison_test(infiltration_loam_sand_1d_u
                           INPUT amanzi_infiltration_loam_sand_1d-u.xml
			   ABSOLUTE_TOLERANCE 1e-7
			   OBSERVATION "aqueous pressure"
			   REFERENCE ${INFILTRATION_LOAM_SAND_1D_SOURCE_DIR}/golden_output/observations.out
			   NPROCS 1)

# add_amanzi_comparison_test(infiltration_loam_sand_1d_s
#                            INPUT amanzi_infiltration_loam_sand_1d-s.xml
# 			     ABSOLUTE_TOLERANCE 1e-8
# 			     OBSERVATION "aqueous pressure"
# 			     REFERENCE observations.out
# 			     NPROCS 1)
