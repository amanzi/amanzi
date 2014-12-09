# -*- mode: cmake -*-

project (BUTLER_STRIP_2D)

set(test_name butler_strip_2d)
set(compfile observations.out)
set(dataset "")
set(step 1)

if (ENABLE_TESTS)
    add_amanzi_test(run_${test_name} KIND AMANZI NPROCS 4 AMANZI_INPUT ${BUTLER_STRIP_2D_SOURCE_DIR}/amanzi_${test_name}.xml)
    add_test(NAME compare_${test_name} 
	     COMMAND diff ${BUTLER_STRIP_2D_BINARY_DIR}/${compfile}
	                  ${BUTLER_STRIP_2D_SOURCE_DIR}/golden_output/${compfile} 
			  #${dataset}/${step}  ${dataset}/${step} 
			  )
    set_tests_properties( compare_${test_name} PROPERTIES DEPENDS run_${test_name} )
    #set_tests_properties( compare_${test_name} PROPERTIES FAIL_REGULAR_EXPRSSION "differences found" )
    set_tests_properties( compare_${test_name} PROPERTIES FAIL_REGULAR_EXPRSSION "differs from" )
    _append_test_label(compare_${test_name} REG)
    _append_test_label(run_${test_name} REG)
endif()
