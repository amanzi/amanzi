# -*- mode: cmake -*-

project (THEIS_ISOTROPIC_1D)

set(test_name theis_isotropic_1d)
set(compfile observation.out)
set(dataset "")
set(step 1)

if (ENABLE_TESTS)
    add_amanzi_test(run_${test_name} KIND AMANZI AMANZI_INPUT ${THEIS_ISOTROPIC_1D_SOURCE_DIR}/amanzi_${test_name}-u.xml NPROCS 4)
    add_test(NAME compare_${test_name} 
	COMMAND ${PYTHON_EXECUTABLE} ${THEIS_ISOTROPIC_1D_SOURCE_DIR}/convergence_test.py ${THEIS_ISOTROPIC_1D_SOURCE_DIR})
	 #COMMAND diff ${THEIS_ISOTROPIC_1D_BINARY_DIR}/${compfile}
	 #                 ${THEIS_ISOTROPIC_1D_SOURCE_DIR}/golden_output/${compfile} 
	 #			  #${dataset}/${step}  ${dataset}/${step} 
	 #			  )
    set_tests_properties( compare_${test_name} PROPERTIES DEPENDS run_${test_name} )
    set_tests_properties( compare_${test_name} PROPERTIES FAIL_REGULAR_EXPRSSION "Comparison Failed" )
    #set_tests_properties( compare_${test_name} PROPERTIES FAIL_REGULAR_EXPRSSION "differs from" )
    _append_test_label(compare_${test_name} REG)
    _append_test_label(run_${test_name} REG)
endif()

