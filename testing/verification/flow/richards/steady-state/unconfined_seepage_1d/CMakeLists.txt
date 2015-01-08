# -*- mode: cmake -*-

project (UNCONFINED_SEEPAGE_1D)

set(test_name unconfined_seepage_1d)
set(compfile observation-isv2.out)
set(dataset "")
set(step 1)

if (ENABLE_TESTS)
    file(COPY porflow4_4.exo DESTINATION .)
    add_amanzi_test(run_${test_name} KIND AMANZI AMANZI_INPUT ${UNCONFINED_SEEPAGE_1D_SOURCE_DIR}/amanzi_${test_name}-isv2.xml)
    add_test(NAME compare_${test_name} 
	     COMMAND diff 
	             ${UNCONFINED_SEEPAGE_1D_BINARY_DIR}/${compfile} 
	             ${UNCONFINED_SEEPAGE_1D_SOURCE_DIR}/golden_output/${compfile} 
		     #${dataset}/${step}  ${dataset}/${step} 
		     )
    set_tests_properties( compare_${test_name} PROPERTIES DEPENDS run_${test_name} )
    #set_tests_properties( compare_${test_name} PROPERTIES FAIL_REGULAR_EXPRSSION "differences found" )
    set_tests_properties( compare_${test_name} PROPERTIES FAIL_REGULAR_EXPRSSION "differs from" )
    _append_test_label(compare_${test_name} REG)
    _append_test_label(run_${test_name} REG)
endif()
