# -*- mode: cmake -*-

project (LINEAR_MATERIALS_SERIAL_1D)

add_amanzi_comparison_test(linear_materials_serial_1d
                           INPUT amanzi_linear_materials_serial_1d-u.xml
			   ABSOLUTE_TOLERANCE 1e-7
			   OBSERVATION "aqueous pressure"
			   REFERENCE ${LINEAR_MATERIALS_SERIAL_1D_SOURCE_DIR}/golden_output/observation.out
			   NPROCS 1)
