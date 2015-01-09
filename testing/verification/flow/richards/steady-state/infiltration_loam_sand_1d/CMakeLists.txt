# -*- mode: cmake -*-

project (INFILTRATION_LOAM_SAND_1D)

set(test_name infiltration_loam_sand_1d)
set(h5file case_2a_plot_data.h5)
set(dataset "/hydraulic_head.cell.0")
set(step 384)
set(gstep 365)

if (ENABLE_TESTS)
    add_amanzi_test(run_${test_name} KIND AMANZI AMANZI_INPUT ${INFILTRATION_LOAM_SAND_1D_SOURCE_DIR}/amanzi_${test_name}-isv2.xml)
    add_test(NAME compare_${test_name} 
	COMMAND ${PYTHON_EXECUTABLE} ${INFILTRATION_LOAM_SAND_1D_SOURCE_DIR}/convergence_test.py ${INFILTRATION_LOAM_SAND_1D_SOURCE_DIR})
	     #COMMAND ${HDF5_H5DIFF_BINARY} ${INFILTRATION_LOAM_SAND_1D_BINARY_DIR}/${h5file} 
	     #                                  ${INFILTRATION_LOAM_SAND_1D_SOURCE_DIR}/golden_output/${h5file} 
	     #					   ${dataset}/${step}  ${dataset}/${gstep} )
    set_tests_properties( compare_${test_name} PROPERTIES DEPENDS run_${test_name} )
    set_tests_properties( compare_${test_name} PROPERTIES FAIL_REGULAR_EXPRSSION "Comparison Failed" )
    #set_tests_properties( compare_${test_name} PROPERTIES FAIL_REGULAR_EXPRSSION "differences found" )
    _append_test_label(compare_${test_name} REG)
    _append_test_label(run_${test_name} REG)
endif()
