# -*- mode: cmake -*-

project (UNCONFINED_RECHARGE_1D)

set(test_name unconfined_recharge_1d)
set(h5file /steady-flow_data.h5)
set(dataset "/hydraulic_head.cell.0")
set(step 3786)

if (ENABLE_TESTS)
    add_amanzi_test(run_${test_name} KIND AMANZI AMANZI_INPUT ${UNCONFINED_RECHARGE_1D_SOURCE_DIR}/amanzi_${test_name}-isv2.xml)
    add_test(NAME compare_${test_name} 
	     COMMAND ${HDF5_H5DIFF_BINARY} ${UNCONFINED_RECHARGE_1D_BINARY_DIR}/${h5file} 
	                                   ${UNCONFINED_RECHARGE_1D_SOURCE_DIR}/golden_output/${h5file} 
					   ${dataset}/${step}  ${dataset}/${step} )
    set_tests_properties( compare_${test_name} PROPERTIES DEPENDS run_${test_name} )
    set_tests_properties( compare_${test_name} PROPERTIES FAIL_REGULAR_EXPRSSION "differences found" )
    _append_test_label(compare_${test_name} REG)
    _append_test_label(run_${test_name} REG)
endif()
