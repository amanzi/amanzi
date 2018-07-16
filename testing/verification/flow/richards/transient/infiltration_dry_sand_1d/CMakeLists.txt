# -*- mode: cmake -*-

project (INFILTRATION_DRY_SAND_1d)

set(test_name infiltration_dry_sand_1d)
set(h5file plot_data.h5)
set(dataset "/hydraulic_head.cell.0")
set(step 1816)

if (ENABLE_TESTS)
    add_amanzi_test(run_${test_name} KIND AMANZI AMANZI_INPUT ${INFILTRATION_DRY_SAND_1d_SOURCE_DIR}/amanzi_${test_name}-u.xml)
    add_test(NAME compare_${test_name} 
	COMMAND ${PYTHON_EXECUTABLE} ${INFILTRATION_DRY_SAND_1d_SOURCE_DIR}/convergence_test.py ${INFILTRATION_DRY_SAND_1d_SOURCE_DIR})
	 #COMMAND ${HDF5_H5DIFF_BINARY} ${INFILTRATION_DRY_SAND_1d_BINARY_DIR}/${h5file}
	 #                                  ${INFILTRATION_DRY_SAND_1d_SOURCE_DIR}/golden_output/${h5file} 
	 #					   ${dataset}/${step}  ${dataset}/${step} )
    set_tests_properties( compare_${test_name} PROPERTIES DEPENDS run_${test_name} )
    set_tests_properties( compare_${test_name} PROPERTIES FAIL_REGULAR_EXPRSSION "Comparison Failed" )
    #set_tests_properties( compare_${test_name} PROPERTIES FAIL_REGULAR_EXPRSSION "differences found" )
    _append_test_label(compare_${test_name} REG)
    _append_test_label(run_${test_name} REG)
endif()
