# -*- mode: cmake -*-

project (BOUNDEDDOMAIN_2D)

set(test_name boundedDomain_2d)
set(compfile observations.out)
set(dataset "")
set(step 1)

if (ENABLE_TESTS)
    add_amanzi_test(run_${test_name} KIND AMANZI AMANZI_INPUT ${BOUNDEDDOMAIN_2D_SOURCE_DIR}/amanzi_${test_name}.xml NPROC 4)
    add_test(NAME compare_${test_name} 
	     COMMAND diff ${BOUNDEDDOMAIN_2D_BINARY_DIR}/${compfile}
	                  ${BOUNDEDDOMAIN_2D_SOURCE_DIR}/golden_output/${compfile} 
			  #${dataset}/${step}  ${dataset}/${step} 
			  )
    set_tests_properties( compare_${test_name} PROPERTIES DEPENDS run_${test_name} )
    #set_tests_properties( compare_${test_name} PROPERTIES FAIL_REGULAR_EXPRSSION "differences found" )
    set_tests_properties( compare_${test_name} PROPERTIES FAIL_REGULAR_EXPRSSION "differs from" )
    _append_test_label(compare_${test_name} REG)
    _append_test_label(run_${test_name} REG)
endif()
