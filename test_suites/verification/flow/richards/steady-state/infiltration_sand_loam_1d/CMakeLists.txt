# -*- mode: cmake -*-

project (INFILTRATION_SAND_LOAM_1D)

add_amanzi_comparison_test(infiltration_sand_loam_1d_u
                           INPUT amanzi_infiltration_sand_loam_1d-u.xml
                           ABSOLUTE_TOLERANCE 1e-8
                           RELATIVE_TOLERANCE 1e-8
			   OBSERVATION "aqueous pressure"
			   REFERENCE ${INFILTRATION_SAND_LOAM_1D_SOURCE_DIR}/golden_output/observations.out
			   NPROCS 1)

# add_amanzi_comparison_test(infiltration_sand_loam_1d_s
#                            INPUT amanzi_infiltration_sand_loam_1d-s.xml
# 			     ABSOLUTE_TOLERANCE 1e-8
# 			     OBSERVATION "aqueous pressure"
# 			     REFERENCE observations.out
# 			     NPROCS 1)
